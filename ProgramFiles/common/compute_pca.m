function [Q,lambdas,Xmean,lambdas_all_sum,my_rel_err] = compute_pca(X,m_max,rel_err)

T = size(X,2);
Xmean = sum(X,2)/T; % compute mean
X = X - kron(ones(1,T),Xmean); % substract mean
mycov = X*X'; % compute covariance matrix
[Q,D] = eigs(mycov,m_max); % compute eigenvalues

lambdas = diag(D);
lambdas_all_sum = trace(mycov);

[lambdas,idx] = sort(lambdas,'descend');
Q = Q(:,idx);

if rel_err > 0
    for new_m = 1:m_max
        my_rel_err = (lambdas_all_sum - sum(lambdas(1:new_m)))/lambdas_all_sum;
        if my_rel_err < rel_err
            break;
        end
    end
    lambdas = lambdas(1:new_m);
    Q = Q(:,1:new_m);
else
    my_rel_err = (lambdas_all_sum - sum(lambdas(1:m_max)))/lambdas_all_sum;
end

end