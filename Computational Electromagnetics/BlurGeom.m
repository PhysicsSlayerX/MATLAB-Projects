% DEFINE GRID
Sx = 1; %physical size along x
Sy = 1; %physical size along y
Nx = 21; %number of cells along x
Ny = 21; %number of cells along y
% GRID ARRAYS
dx = Sx/Nx;
xa = [0:Nx-1]*dx;
xa = xa - mean(xa);
dy = Sy/Ny;
ya = [0:Ny-1]*dy;
ya = ya - mean(ya);
[Y,X] = meshgrid(ya,xa);
% CREATE A CROSS
ER = abs(X)<=0.075 | abs(Y)<=0.075;
% CREATE BLUR FUNCTION
B = exp(-(X.^2 + Y.^2)/0.1^2);
% PERFORM BLUR OPERATION
ER = fft2(ER).*fft2(B)/sum(B(:));
ER = ifftshift(real(ifft2(ER)));
% PERFORM THRESHOLD OPERATION
ER = ER > 0.4;