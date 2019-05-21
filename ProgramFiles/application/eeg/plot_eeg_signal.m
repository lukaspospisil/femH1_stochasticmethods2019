function plot_eeg_signal( X, X_rec, R_to_plot )

T = size(X,2);

figure
for r_idx = 1:length(R_to_plot)
    r = R_to_plot(r_idx);
    
    subplot(length(R_to_plot),1,r_idx)
%    title([num2str(r) '. electrode'])

    hold on
    plot(1:T,X(r,:),'b')
    plot(1:T,X_rec(r,:),'r')
    axis([1,T,-200,200])

%    xlabel('$t$','Interpreter','latex')
    ylabel(['$x_{' num2str(r) '}(t)$'],'Interpreter','latex')

    % smart axis (do not repeat x-axis for all subfigures)
    if r_idx == 1
        legend('original', 'recovered')
    end
    if r_idx == length(R_to_plot)
        xlabel('$t$','Interpreter','latex')
    end
    
    hold off
end



end

