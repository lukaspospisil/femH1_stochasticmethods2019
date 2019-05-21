function [ aa ] = norm_fro3( A )

% frobenius norm for 3D tensor
aa = 0;
for i=1:size(A,3)
    aa = aa + norm(A(:,:,i),'fro')^2;
end

aa = sqrt(aa);

end

