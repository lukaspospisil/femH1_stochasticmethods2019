function [ X_orig ] = generate_signal_n5_K3( Trepeat )

% affiliation for one period
aff_one_period = [1*ones(1,100) 2*ones(1,80) 1*ones(1,50) 3*ones(1,70) 2*ones(1,90) 1*ones(1,50) 3*ones(1,60)];
K = max(aff_one_period);

T = length(aff_one_period);
Gamma = zeros(K,T);
for k=1:K
    Gamma(k,aff_one_period == k) = 1;
end

Theta = [1, 0, 2, -1, 0; ...
         2, 0, -1, 0, 1; ...
         -1, 0, 2, 0, 1.5]';
     
X_orig = kron(ones(1,Trepeat),Theta*Gamma);

end

