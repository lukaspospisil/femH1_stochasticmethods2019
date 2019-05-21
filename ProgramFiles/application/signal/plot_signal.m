function plot_signal( X_true, X, X_rec )

n = max([size(X_true,1),size(X,1),size(X_rec,1)]);

% if something is [], then do not plot it
mylegend = {}; % here store the legend of figure

figure
for i=1:n
    subplot(n,1,i);
    hold on
    
    if ~isempty(X)
        plot(1:size(X,2),X(i,:),'b','Color',[0.7,0.7,0.7],'LineWidth',1.0);
        mylegend{numel(mylegend)+1} = 'data';
    end
    
    if ~isempty(X_true)
        plot(1:size(X_true,2),X_true(i,:),'g','Color',[0,0.6,0],'LineWidth',2.0);
        mylegend{numel(mylegend)+1} = 'true';
    end
    
    if ~isempty(X_rec)
        plot(1:size(X_rec,2),X_rec(i,:),'r','LineWidth',2.0);
        mylegend{numel(mylegend)+1} = 'reconstructed';
    end
    
    xlabel('$t$','Interpreter','latex');
    ylabel(['$x_{' num2str(i) '}(t)$'],'Interpreter','latex');
    
    if i==1
        legend(mylegend);
    end
    
    hold off
end

end

