function [abserr_min, idx_min] = plot_abserr( epssqrs,abserrs )

figure
hold on
plot(epssqrs,abserrs,'ro-','linewidth',1.5)
xlabel('$\epsilon^2$','interpreter','latex')
ylabel('absolute error')
set(gca,'xscale','log')

% find the minimum
[~,idxs] = min(abserrs);
idx_min = idxs(end);
abserr_min = abserrs(idx_min);

% add text
% but maybe only the first and the last epssqr
plot(epssqrs(idx_min),abserr_min,'b*')

text(epssqrs(idx_min),abserr_min,['$\epsilon = ' num2str(epssqrs(idx_min)) '$'],'interpreter','latex')

hold off

end

