function [ out ] = signalclustering( in )

% in.X in format [n,T]
[n,T] = size(in.X);

% if S is given, then annealing doesn't make any sence - Gamma problem is
% convex
if ~isempty(in.options.S_given)
    in.options.nanneal = 1;
end

Hreg = get_laplacian_grid1d(T); % time regularization matrix

out.L = Inf;

for idx_ann = 1:in.options.nanneal
    if strcmp(in.options.type,'FEMH1_quadprog')
        [ S, Gamma, L, Llin, Lquad ] = femH1_quadprog(in.X, Hreg, in.options);
    end
    if strcmp(in.options.type,'FEMH1')
        [ S, Gamma, L, Llin, Lquad ] = femH1(in.X, Hreg, in.options);
    end
    if strcmp(in.options.type,'FEMH1_gpu')
        [ S, Gamma, L, Llin, Lquad ] = femH1_gpu(in.X, Hreg, in.options);
    end

    disp(['ann = ' num2str(idx_ann) ', L = ' num2str(L)])
    
    if L < out.L
        out.L = L;
        out.S = S;
        out.Gamma = Gamma;
        out.Llin = Llin;
        out.Lquad = Lquad;
    end
    
    in.options.Gamma0 = []; % next annealing will have different init approx
end
    
end
