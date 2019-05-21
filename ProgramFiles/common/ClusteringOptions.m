classdef ClusteringOptions
    
    properties
        type            % 'FEMH1'/'FEMH1_quadprog'/'FEMH1_gpu'
        
        K               % number of clusters
        epssqr          % regularization parameter
        maxit           % max number of outer (subspace) algorithm (only if S is not given)
        eps             % stopping criteria of outer (subspace) algorithm (only if S is not given)
        qp_maxit        % max number of QP solver iterations
        qp_eps          % stopping criteria of QP solver
        nanneal         % number of annealing steps
        
        S_given         % given parameters of clusters
        Gamma0          % given initial approximation of subspace algorithm (if empty, then algorithm generates random)
        
        gpudev          % on which GPU to compute, result of gpuDevice() (only for *_gpu variant)
        
        dispdebug       % plot info about progress or be quite (true/false)
        
        pca_m           % PCA parameter
    end
    
    methods
        function obj = ClusteringOptions()
            % constructor - set default values
            obj.type = 'FEMH1';
            
            obj.K = 2;
            obj.epssqr = 1e-2;
            obj.maxit = 100;
            obj.eps = 1e-4;
            obj.qp_maxit = 1e2;
            obj.qp_eps = 1e-4;
            obj.nanneal = 5;
            
            obj.S_given = [];
            obj.Gamma0 = [];
            
            obj.gpudev = [];
            
            opt.dispdebug = true;
            
            opt.pca_m = 5;
        end
        
    end
end    
