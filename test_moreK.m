% compute results for several K=1,...

clear all

addpath(genpath('ProgramFiles')) % add all files in 'ProgramFiles' directory

%% define sample signal
if true
    % load data from file
    X = load('data.mat');
    X = X.X; % sorry for this
    X_true = []; % no, we don't have true signal
else
    % generate new data
    addpath('data/signal') % benchmark generator
    X_true = generate_signal_n1_K3(1);
%    X_true = generate_signal_n5_K3(1);
    sigma = 1.5; % noise parameter
    X = X_true + sigma*randn(size(X_true)); % add noise
end

disp(['Signal: [n,T] = ' num2str(size(X,1)) ' x ' num2str(size(X,2)) ' = ' num2str(prod(numel(X)))]);

% set options of clustering algorithm
in.X = X;
in.options = ClusteringOptions(); % set default options
in.options.type = 'FEMH1'; % FEMH1_quadprog / FEMH1 / FEMH1_gpu
in.options.S_given = [];%[]; % given parameters of clusters (of size n,K), if not given then S_given = []
in.options.qp_eps = 1e-6; % the precision of QP solver
in.options.qp_maxit = 1e2; % max number of iterations of QP solver
in.options.eps = 1e-4; % the precision of subspace algorithm
in.options.maxit = 20; % max number of subspace algorithm
in.options.dispdebug = false; % display some info about progress
in.options.nanneal = 1; % number of annealing steps

% prepare gpu
if strcmp(in.options.type, 'FEMH1_gpu')
    in.options.gpudev = gpuDevice();
    reset(in.options.gpudev);
end

% solve problem for these values of regularization parameter and numbers of
% clusters
Ks = [1,2,3,4,5,6];
%epssqrs = 10.^[-7:0.5:3];
epssqrs = 10.^[-7:1:3];

% here store the results
Ls = zeros(length(Ks),length(epssqrs));
Llins = zeros(length(Ks),length(epssqrs));
Lquads = zeros(length(Ks),length(epssqrs));
abserrs = zeros(length(Ks),length(epssqrs)); % only if X_true is known

% for every K solve the problem
for idx_Ks = 1:length(Ks)
    in.options.K = Ks(idx_Ks); % set number of clusters
    disp(['K = ' num2str(in.options.K) ' (' num2str(idx_Ks) ' of ' num2str(length(Ks)) ')']);

    % for every regularization parameter solve the problem
    for idx_epssqrs = 1:length(epssqrs)
        in.options.epssqr = epssqrs(idx_epssqrs); % regularization parameter
        disp([' - epssqr = ' num2str(in.options.epssqr) ' (' num2str(idx_epssqrs) ' of ' num2str(length(epssqrs)) ')'])

        % solve the problem using our clustering method
        [out] = signalclustering(in); % this is it!

        % store results
        Ls(idx_Ks,idx_epssqrs) = out.L;
        Llins(idx_Ks,idx_epssqrs) = out.Llin;
        Lquads(idx_Ks,idx_epssqrs) = out.Lquad;
        
        % I would like to compute also absolute error of reconstruction
        % but in this case, I will need X_true
        if ~isempty(X_true)
            % compute reconstructed signal
            [ X_rec ] = xrec_signal( out.S, out.Gamma );

            abserrs(idx_Ks,idx_epssqrs) = norm(X_rec - X_true);
        end
    
        % idea - use solution from previous epssqr as init approximation for
        % next epssqr
        in.options.Gamma0 = out.Gamma;
    end
    
    in.options.Gamma0 = []; % new gamma has different dimension
end

% plot results
test_moreK_plot

