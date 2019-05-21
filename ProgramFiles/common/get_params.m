function [ params ] = get_params( Gamma,S )

% compute number of parameters to store Gamma
[K,T] = size(Gamma);

% the eigenvectors of Laplace matrix
k = 0:T-1;
V = normc(cos(pi*kron((1:T)',k)/T - pi*k/(2*T)));

rel_err = 1e-6; % we want to have relative error of reconstruction smaller than this value

params_Gamma = 0;
for idx_Ks = 1:K
   m = V'*Gamma(idx_Ks,:)';
   
   % motivated by PCA
   m_sum = sum(m);
   for new_m = 1:T
        my_rel_err = (m_sum - sum(m(1:new_m)))/m_sum;
    
        if my_rel_err < rel_err
            break;
        end
   end
   params_Gamma = params_Gamma + new_m;
end
       
% compute number of parameters to store S
n = size(S,1);
params_S = K*n;

params = params_Gamma + params_S;

end

