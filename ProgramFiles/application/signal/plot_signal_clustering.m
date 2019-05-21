function plot_signal_clustering( Gamma, S )

[K,T] = size(Gamma);

figure
for k=1:K
    subplot(K,1,k);
    hold on
    if ~isempty(S)
        title(['$S_' num2str(k) ' = ' vec_to_nice_str(S(:,k)) '$'], 'Interpreter', 'latex');
    else
        title(['$S_' num2str(k) '$'], 'Interpreter', 'latex');
    end
    plot(1:T,Gamma(k,:),'b')
    
    xlabel('$t$','Interpreter','latex')
    ylabel(['$\gamma_{' num2str(k) '}(t)$'],'Interpreter','latex')
    
    axis([1 T -0.2 1.2])
    hold off
end

end

function mystr = vec_to_nice_str(v)
% transform vector v = [1,2,3] to string '[1,2,3]'
% transform e-7 to 10^{-7} to have it completely awesome
mystr = '[';
for i = 1:length(v)
    numstr = num2str(v(i));
    
    if contains(numstr, 'e') % ~isempty(strfind(numstr, 'e'))
        numstr = strrep(numstr,'e','\cdot 10^{');
        numstr = [numstr '}'];
    end
    
    mystr = [mystr numstr];
    if i < length(v)
        mystr = [mystr ','];
    end
end
mystr = [mystr ']'];
end
