function [ S, Gamma, L, Llin, Lquad ] = femH1_gpu(X,Hreg,options )

% provide the signal in format [n,T]
[n,T] = size(X);
K = options.K;

% check sizes of input variables
if size(Hreg,1) ~= T
    error('! size of regularization matrix does not match with the size of data')
end
if ~isempty(options.S_given)
    if size(options.S_given,1) ~= n
        error('! n of data and size of S_given do not match')
    end
    if size(options.S_given,2) ~= K
        error('! given K and size of S_given do not match')
    end
end
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

% generate initial random feasible Gamma [K,T] or use provided one
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
S = zeros(n,K); 

% initial object function value
L = Inf;

% move everything to GPU
L = gpuArray(L);
Gamma = gpuArray(Gamma);
gamma_vec = gpuArray(gamma_vec);
g = gpuArray(g);
HG = gpuArray(HG);
S = gpuArray(S);

it = 0; % iteration counter
while it < maxit % practical stopping criteria is present after computing new L (see "break")
    
    % compute S
    if isempty(options.S_given)
        for k=1:K
            sum_gammak = sum(Gamma(k,:));
            if sum_gammak ~= 0 % maybe gamma_k = 0 ? (i.e. this cluster is empty)
                S(:,k) = sum(bsxfun(@times,X,Gamma(k,:)),2)/sum_gammak;
            else
                S(:,k) = zeros(n,1);
            end
        end
    else
        S = options.S_given;
    end
    
    % compute new Gamma
    % prepare new linear term in QP,
    % i.e. compute new residuum based on new S
    g = sum((kron(ones(1,K),X) - kron(S,ones(1,T))).^2,1)';

    % rescale linear term
    g = g/(T*n);

    % solve QP problem
%    [gamma_vec,itQP] = spgqp((2*options.epssqr)*HG, -g, gamma_vec, K, options.epssqr*normH, options.qp_eps, options.qp_maxit);
    [gamma_vec,itQP] = spgqp_gpu((2*options.epssqr)*HG, -g, gamma_vec, K, options.epssqr*normH, options.qp_eps, options.qp_maxit, options.gpudev);

    Gamma = reshape(gamma_vec,T,K)';

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
    Lit(it) = L; % for postprocessing
end

% move results back grom GPU
wait(options.gpudev);
S = gather(S);
Gamma = gather(Gamma);
L = gather(L);
Llin = gather(Llin);
Lquad = gather(Lquad);

end

