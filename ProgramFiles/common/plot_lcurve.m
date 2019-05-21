function [lcurve_max, idx_max] = plot_lcurve( epssqrs,Llins,Lquads )

figure
hold on
title('L-curve')
plot(Llins,Lquads,'bo-','linewidth',1.5)
xlabel('modelling error')
ylabel('regularization')

% add text
% but maybe only the first and the last epssqr
toplot = [1,length(epssqrs)];
for idx=1:length(toplot)
    text(Llins(toplot(idx)),Lquads(toplot(idx)),['$\epsilon = ' num2str(epssqrs(toplot(idx))) '$'],'interpreter','latex');
end


% compute curvature of L-curve
% https://en.wikipedia.org/wiki/Curvature
d_Lquads = diff(Lquads,1);
d_Llins = diff(Llins,1);
dd_Lquads = diff(Lquads,2);
dd_Llins = diff(Llins,2);
d_Lquads = 0.5*(d_Lquads(2:end)+d_Lquads(1:end-1));
d_Llins = 0.5*(d_Llins(2:end)+d_Llins(1:end-1));

if sum(d_Llins.^2 + d_Lquads.^2) > 0
    kappa = abs(d_Llins.*dd_Lquads - dd_Llins.*d_Lquads)./((d_Llins.^2 + d_Lquads.^2).^(3/2));
else
    kappa = zeros(size(dd_Llins));   
end

[~, idx_max] = max([0; kappa(2:end-1); 0]); % find max value of curvature
idx_max = idx_max(end);
lcurve_max = kappa(idx_max);

% plot max
plot(Llins(idx_max),Lquads(idx_max),'r*');
text(Llins(idx_max),Lquads(idx_max),['$\epsilon = ' num2str(epssqrs(idx_max)) '$'],'interpreter','latex');

hold off


end

