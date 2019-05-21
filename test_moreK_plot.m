% run after test_moreK

% compute postprocess through all possible K
K_Ls = zeros(length(Ks),1);
K_Llins = zeros(length(Ks),1);
K_Lquads = zeros(length(Ks),1);
K_abserrs = zeros(length(Ks),1);

% store indexes of "best epssqrs"
K_lcurveopt_idx = zeros(length(Ks),1);
K_abserropt_idx = zeros(length(Ks),1);

for idx_Ks = 1:length(Ks)
    % plot Lcurves
    [lcurve_max, idx_lcurveopt] = plot_lcurve(epssqrs,Llins(idx_Ks,:)',Lquads(idx_Ks,:)');
    title(['L-curve, K = ' num2str(Ks(idx_Ks))])

    K_lcurveopt_idx(idx_Ks) = idx_lcurveopt;
    
    % plot absolute error
    if ~isempty(X_true)
        [abserr_min, idx_abserropt] = plot_abserr(epssqrs,abserrs(idx_Ks,:));
        title(['absolute error curve, K = ' num2str(Ks(idx_Ks))])
        K_abserropt_idx(idx_Ks) = idx_abserropt(end);

        K_abserrs(idx_Ks) = abserr_min;
    end

    idx_opt = K_lcurveopt_idx(idx_Ks);

    K_Ls(idx_Ks) = Ls(idx_Ks,idx_opt);
    K_Llins(idx_Ks) = Llins(idx_Ks,idx_opt);
    K_Lquads(idx_Ks) = Lquads(idx_Ks,idx_opt);
    
end

% plot Ls for different K
figure
hold on
plot(Ks,K_Llins,'bo-','linewidth',1.5)
xlabel('K')
ylabel('modelling error')

