function [ S, Gamma, L, Llin, Lquad ] = femH1pca(X,Hreg,options )

% provide the signal in format [n,T]
[n,T] = size(X);
K = options.K;

%TODO: make as option
m_max = 2;
rel_err = 0;%1e-12;

% check sizes of input variables
if size(Hreg,1) ~= T
    error('! size of regularization matrix does not match with the size of data')
end
% TODO: check options.S_given
if ~isempty(options.Gamma0)
    if size(options.Gamma0,1) ~= K
        error('! K and size of Gamma0 do not match')
    end
    if size(options.Gamma0,2) ~= T
        error('! length of data and Gamma0 do not match')
    end
end

% if someone provided S, then there are not outer iterations
if isempty(options.S_given)
    maxit = options.maxit;
else
    maxit = 1;
end

% generate initial random feasible Gamma [K,T]
if ~isempty(options.Gamma0)
    Gamma = options.Gamma0;
else
    Gamma = rand(K,T);
end

% initial gamma has to be feasible (i.e sum to one for each t and positive)
% this is cheaper than projection
Gamma_sum = sum(Gamma,1);
for k=1:K
    Gamma(k,:) = Gamma(k,:)./Gamma_sum;
end
gamma_vec = reshape(Gamma',T*K,1); % vectorized form of gamma used in QP

% prepare QP objects
HG = kron(speye(K),Hreg);

% in the case of 1D laplace, we have analytical formula for the largest eig
normH = (2 - 2*cos(pi*k/T));
%normH = gershgorin(Hreg); % estimate largest eigenvalue

g = zeros(size(gamma_vec)); % linear term in QP

% here will be stored solution of model parameters - one for each cluster
S.Q = cell(K,1);
S.lambdas = cell(K,1);
S.Xmean = cell(K,1);
S.coeff = cell(K,1);
S.lambdas_all_sum = cell(K,1);

% initial object function value
L = Inf;

it = 0; % iteration counter
while it < maxit % practical stopping criteria is present after computing new L (see "break")
    
    % compute S
    if isempty(options.S_given)
        Gamma2 = round_gamma(Gamma);
        for k=1:K
            if sum(Gamma2(k,:)) > 0
                [S.Q{k},S.lambdas{k},S.Xmean{k},S.lambdas_all_sum{k},S.rel_err{k}] = ...
                    compute_pca(X(:,Gamma2(k,:) == 1),m_max,rel_err);
                
                Xreduced = S.Q{k}'*(X - kron(ones(1,T),S.Xmean{k}));
                Xrec = S.Q{k}*Xreduced + kron(ones(1,T),S.Xmean{k});
                
                g((k-1)*T+1:k*T) = sum((X - Xrec).^2,1);
            else
                g((k-1)*T+1:k*T) = 0;
                
            end
        end
        
    else
        S = options.S_given; % S is provided - there is nothing to solve
    end
    
    % rescale linear term (to avoid huge numbers if T is large)
    g = g/(T*n);
    
    % solve QP problem
    [gamma_vec,itQP] = spgqp((2*options.epssqr)*HG, -g, gamma_vec, K, options.epssqr*normH, options.qp_eps, options.qp_maxit);
    
    Gamma = reshape(gamma_vec,T,K)'; % TODO: it is really necessary to keep both of Gamma and gamma_vec?
    
    % recompute coeffs based on new gamma
    Gamma2 = round_gamma(Gamma);
    for k=1:K
        if and(sum(Gamma2(k,:)) > 0,~isempty(S.Q{k}))
            S.coeff{k} = S.Q{k}'*(X(:,Gamma2(k,:)==1) - kron(ones(1,sum(Gamma2(k,:))),S.Xmean{k}));
        end
    end
    
    % compute new function value
    Lold = L; % store old function value
    Llin = dot(g,gamma_vec);
    Lquad = dot(HG*gamma_vec,gamma_vec);
    L = Llin + options.epssqr*Lquad;
    deltaL = Lold - L; % the progress of objective function, Lold > L (?)
    
    % display some informations, so we know that something is happening
    if options.dispdebug
        disp([num2str(it) '. it: L = ' num2str(L) ', deltaL = ' num2str(deltaL) ', itQP = ' num2str(itQP)])
    end
    
    % stopping (breaking) criteria
    % based on sufficient decrease of objective funtion
    if abs(deltaL) < options.eps
        break; % stop outer "while" cycle
    end
    
    it = it + 1; % increase number of iterations
    Lit(it) = L; % for postprocessing & fun
end


end

