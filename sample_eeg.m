clear all

addpath(genpath('ProgramFiles')) % add all files in 'ProgramFiles' directory

%% load data
[hdr, X] = edfRead('data/eeg/S001R01.edf.txt');
X_true = []; % we don't have exact signal

% create a matrix of EEG time series, the rows correspond
% to 64 measurement electrodes (numbered from 1 to 64),
% the columns correspond to different equidistant time measurements
% For the sake of simplicity and speed, we will first consider only
% the first few measurements from these data
%in.X=X(1:64,1:500);
in.X=X(1:64,:);

% print informations of the data
[R,T]=size(in.X);
disp(['loaded data: R = ' num2str(R) ', T = ' num2str(T)])

% load the coordinates of 64 electrode positions (in variable P,
% first column of P refers to the coordinates in x-direction,
% second column of P refers to the coordinates in y-direction )
load data/eeg/Koordinaten_EEG2

in.P = P; % provide graph information (vertices)

%% solve clustering problem
% set options of clustering algorithm
in.options = ClusteringOptions(); % set default options
in.options.type = 'FEMH1PCA_gpu'; % FEMH1_quadprog / FEMH1 / FEMH1_gpu
in.options.K = 3; % set number of clusters
in.options.epssqr = 1e1; % regularization parameter
in.options.S_given = []; % given parameters of clusters (of size n,K), if not given then S_given = []
in.options.qp_eps = 1e-12; % the precision of QP solver
in.options.qp_maxit = 2e2; % max number of iterations of QP solver
in.options.eps = 1e-6; % the precision of subspace algorithm
in.options.maxit = 30; % max number of subspace algorithm
in.options.dispdebug = true; % display some info about progress or be quite (true/false)
in.options.nanneal = 1;

% prepare gpu
if strcmp(in.options.type, 'FEMH1PCA_gpu')
    in.options.gpudev = gpuDevice();
    reset(in.options.gpudev);
    
    addpath('CudaFiles')
end

% solve the problem using our clustering method
[out] = eegclustering(in);

% it is k-means like analysis
if or(strcmp(in.options.type,'FEMH1'),strcmp(in.options.type,'FEMH1_gpu'))
    % compute reconstructed signal
    [ X_rec ] = xrec_signal( out.S, out.Gamma );

    % plot reconstruction
    if true
        R_to_plot = [1,2,3,4,5]; % which electrodes to plot

        if isempty(X_true)
            % do not plot true signal - we don't have it
            % compare just data and reconstruction
            plot_signal( [], X(R_to_plot,:), X_rec(R_to_plot,:) );
        else
            % compare true signal, data, and reconstruction
            plot_signal( X_true(R_to_plot,:), X(R_to_plot,:), X_rec(R_to_plot,:) );
        end
    end
    
    % plot Gamma - clustering affiliation
    if true
        plot_signal_clustering( out.Gamma, [] );
    end
end

if strcmp(in.options.type,'FEMH1PCA')
    % compute reconstructed signal
    [ X_rec ] = xrec_eeg( out.S, out.Gamma );
    
    if true
        R_to_plot = [1,2,3,4,5]; % which electrodes to plot
        plot_eeg_signal( in.X, X_rec, R_to_plot );
        % plot Gamma - clustering affiliation
        plot_signal_clustering( out.Gamma, [] );
    end
    
    if false
        T_to_plot = 1:10;
        plot_eeg_movie( in.P, M, X_rec, T_to_plot );
    end
end