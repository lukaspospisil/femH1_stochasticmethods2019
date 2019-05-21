% testing the scaling of QP solvers fro nice picture in presentation

clear all

addpath(genpath('ProgramFiles')) % add all files in 'ProgramFiles' directory
addpath('data/signal') % benchmark generator

in.options = ClusteringOptions();
in.options.K = 3; % set number of clusters
in.options.epssqr = 5e-2; % regularization parameter
in.options.S_given = [1,2,3];%[]; % given parameters of clusters, if not given then S_given = []
in.options.qp_eps = 1e-6; % the precision of QP solver
in.options.qp_maxit = 1e2; % max number of iterations of QP solver
in.options.eps = 1e-4; % the precision of subspace algorithm
in.options.maxit = 20; % max number of subspace algorithm

in.options.gpudev = gpuDevice();

sigma = 0.5; % noise parameter

Trepeats = [1:5,10:5:50,50:10:200,250:50:1000,2000:1000:10000];%10:10:200;
Ts = zeros(length(Trepeats),1);

methods = {'FEMH1_quadprog','FEMH1','FEMH1_gpu'}; % test these algorithms
mylegend = {'FEMH1 quadprog','FEMH1 SPG-QP','FEMH1 SPG-QP GPU'}; % what to show in fig legend

Tmax = [50,5000,Inf]; % max Trepeat for each method (to not wait ages)

times = Inf*ones(length(Trepeats),numel(methods)); % computational time for each method
Ls = Inf*ones(size(times));

for idx_Ts = 1:length(Trepeats)

    %% generate sample signal
    X_orig = generate_signal_n1_K3(Trepeats(idx_Ts));

    disp(['Signal: [n,T] = ' num2str(size(X_orig,1)) ' x ' num2str(size(X_orig,2)) ' = ' num2str(prod(numel(X_orig)))]);
    Ts(idx_Ts) = size(X_orig,2);
    
    % add noise
    in.X = X_orig + sigma*randn(size(X_orig));

    % generate initial Gamma used for all algorithms
    in.Gamma0 = rand(in.options.K,size(X_orig,2));
    
    % solve using different methods
    for idx_method=1:numel(methods)
        in.options.type = methods{idx_method};
    
        if Trepeats(idx_Ts) <= Tmax(idx_method)
            tic_clustering = tic;
            [out] = signalclustering(in);
            times(idx_Ts,idx_method) = toc(tic_clustering);
            
            Ls(idx_Ts,idx_method) = out.L;
        end
        
    end
end

figure
hold on
plot(Ts,times,'o-')
xlabel('T')
ylabel('time [s]')
legend(mylegend)
hold off

figure
hold on
plot(Ts,Ls,'o-')
xlabel('T')
ylabel('L')
legend(mylegend)
hold off