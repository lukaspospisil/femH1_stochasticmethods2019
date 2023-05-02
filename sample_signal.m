clear all

addpath(genpath('ProgramFiles')) % add all files in 'ProgramFiles' directory

%% define sample signal
if false
    % load data from file
    X = load('data.mat');
    X = X.X; % sorry for this
    X_true = []; % no, we don't have true signal
else
    % generate new data
    addpath('data/signal') % benchmark generator
    X_true = generate_signal_n1_K3(10);
%    X_true = generate_signal_n5_K3(1);
    sigma = 0.5; % noise parameter
    X = X_true + sigma*randn(size(X_true)); % add noise
%    save('data.mat','X'); % save the data, maybe we will use it later?
end

disp(['Signal: [n,T] = ' num2str(size(X,1)) ' x ' num2str(size(X,2)) ' = ' num2str(prod(numel(X)))]);

%% solve clustering problem
% set options of clustering algorithm
in.X = X;
in.options = ClusteringOptions(); % set default options
in.options.type = 'FEMH1'; % FEMH1_quadprog / FEMH1 / FEMH1_gpu
in.options.K = 3; % set number of clusters
in.options.epssqr = 1e-2; % regularization parameter
in.options.S_given = [1,2,3];%[]; % given parameters of clusters (of size n,K), if not given then S_given = []
in.options.qp_eps = 1e-6; % the precision of QP solver
in.options.qp_maxit = 1e2; % max number of iterations of QP solver
in.options.eps = 1e-4; % the precision of subspace algorithm
in.options.maxit = 30; % max number of subspace algorithm
in.options.dispdebug = true; % display some info about progress or be quite (true/false)

% prepare gpu
if strcmp(in.options.type, 'FEMH1_gpu')
    in.options.gpudev = gpuDevice();
    reset(in.options.gpudev);
    
    addpath('CudaFiles')
end

% solve the problem using our clustering method
tic_clustering = tic;
[out] = signalclustering(in); % this is it!
time_clustering = toc(tic_clustering);
disp(['Problem solved in ' num2str(time_clustering) 's'])

% compute reconstructed signal
[ X_rec ] = xrec_signal( out.S, out.Gamma );

% plot reconstruction
if true
    if isempty(X_true)
        % do not plot true signal - we don't have it
        % compare just data and reconstruction
        plot_signal( [], X, X_rec );
    else
        % compare true signal, data, and reconstruction
        plot_signal( X_true, X, X_rec );
    end
end

% plot Gamma - clustering affiliation
if true
    plot_signal_clustering( out.Gamma, out.S );
end

