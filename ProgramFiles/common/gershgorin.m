function [ normA ] = gershgorin( A )

% estimation of the largest matrix eigenvalue using Gershgorin circle
% theorem - vectorized variant, works also on GPU
normA_row = diag(A) - abs(diag(A)) + sum(abs(A),2);
normA = max(normA_row);

end

