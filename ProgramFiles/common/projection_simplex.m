function PX = projection_simplex(X,K)

% input X is a vector of size KT
T = length(X)/K;

XX = reshape(X,T,K)';

PXX = projection_simplexes(XX);

PX = reshape(PXX',K*T,1);

end