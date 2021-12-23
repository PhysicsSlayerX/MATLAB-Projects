% figure
% P = peaks(20);
% X = repmat(P, [5,10]);
% imagesc(X);
% colorbar;
% 
% Y = fft2(X).*exp(X);
% imagesc(abs(fftshift(Y)));
% colorbar;
% 
% 
% Z = ifftshift(real(ifft2(Y)));
% imagesc(Z);
% colorbar;

rng(100);
A = randi(10,3);
A_P = A(:)
%B = reshape(A_P, 3, 3)

A + 3