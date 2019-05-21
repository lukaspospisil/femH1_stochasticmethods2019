function [ A ] = get_laplacian_grid1d( T )

% get 3-diag Laplace matrix
if T > 1
    ondiag = [ 1 2*ones(1,T-2) 1];
    offdiag = -ones(1,T);
    A = spdiags([offdiag' ondiag' offdiag'],-1:1,T,T);
else
    A = 0;
end

end

