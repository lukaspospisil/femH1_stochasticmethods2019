function [ Xrec ] = xrec_eeg( S, Gamma )

[K,T] = size(Gamma);
R = size(S.Q{1},1);

% compute reconstructed data
Xrec = zeros(R,T);
Gamma = round_gamma(Gamma);
for k=1:K
    if sum(Gamma(k,:)) > 0
        Xrec(:,Gamma(k,:)==1) = S.Q{k}*S.coeff{k} + kron(ones(1,sum(Gamma(k,:))),S.Xmean{k});
    end
end

end

