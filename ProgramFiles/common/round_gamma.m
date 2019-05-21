function [Gamma] = round_gamma(Gamma_in)

[~,max_idx] = max(Gamma_in,[],1);
Gamma = zeros(size(Gamma_in));
for k=1:size(Gamma_in,1)
    Gamma(k,max_idx==k) = 1;
end

end