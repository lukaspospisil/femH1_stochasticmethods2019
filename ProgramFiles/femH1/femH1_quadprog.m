function [ S, Gamma, L, Llin, Lquad ] = femH1_quadprog(X,Hreg,options )

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


% if someone provides S or Gamma, then there are not outer iterations
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
Gamma_sum = sum(Gamma,1);
for k=1:K
    Gamma(k,:) = Gamma(k,:)./Gamma_sum;
end
gamma_vec = reshape(Gamma',T*K,1); % vectorized form of gamma used in QP

% prepare QP objects for quadprog
HG = kron(speye(K),Hreg); % Hessian matrix (constant during iterations)
g = zeros(size(gamma_vec)); % this will be linear term in QP (changing in every outer it)

% create equality constraints (constant in outer it)
B = kron(ones(1,K),speye(T));
c = ones(T,1);

% lower bound
l = zeros(K*T,1);

% settings of algorithm (quadprog)
options_quadprog = optimoptions('quadprog','Display','none', ... % set 'Display' to 'iter' if you want to see the progress of quadprog
        'Algorithm','interior-point-convex', ...
        'OptimalityTolerance',options.qp_eps, ...
        'MaxIterations',options.qp_maxit);

% here will be stored solution of model parameters - one for each cluster
S = zeros(n,K); 

% initial object function value
L = Inf;

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

    % solve QP problem using quadprog
    [gamma_vec, ~,~, output_quadprog] = quadprog((2*options.epssqr)*HG,g,[],[],B,c,l,[],gamma_vec,options_quadprog);
    itQP = output_quadprog.iterations;
    
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
        disp(['     err_eq = ' num2str(norm(B*gamma_vec-c)) ', err_ineq = ' num2str(norm(min(gamma_vec,0))) ', err_ineq2 = ' num2str(norm(max(gamma_vec-1,0)))])
    end
    
    % stopping (breaking) criteria
    % based on sufficient decrease of objective funtion
    if abs(deltaL) < options.eps
        break; % stop outer "while" cycle
    end
    
    it = it + 1; % increase number of iterations
    Lit(it) = L; % for postprocessing
end


end

