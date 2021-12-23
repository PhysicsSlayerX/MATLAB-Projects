% demo_TO.m
% Implementation of Transformation Optics

% INITIALIZE MATLAB
close all;
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBJECTS
w = 0.65;            % dimention of square cloak
t = 0.4;            % length of triangle side

% GRID
Sx = 1;
Sy = Sx;
% Nx = 100;
% Ny = round(Nx*Sy/Sx);

Nx = 840;
Ny = 700;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1 -- BUILD CLOAK AND OBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID
dx = Sx/Nx;
dy = Sy/Ny;
xa = [0.5:Nx-0.5]*dx;       xa = xa - mean(xa);
ya = [0.5:Ny-0.5]*dy;       ya = ya - mean(ya);

% BUILD TRIANGLE CLOAK
h = w*sqrt(3)/2;
ny = round(w/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;
CLK = zeros(Nx, Ny);

for ny = ny1 : ny2
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f * w / dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    CLK(nx1:nx2,ny) = 1;
end

% BUILD TRIANGLE
h = t*sqrt(3)/2;
ny = round(t/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;
OBJ = zeros(Nx, Ny);

for ny = ny1 : ny2
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f * t / dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    OBJ(nx1:nx2,ny) = 1;
end

% SHOW CLOAK AND OBJECT
subplot(1,2,1);
imagesc(xa, ya, CLK.');
colorbar;
axis equal tight;
title('CLOAK')

subplot(1,2,2);
imagesc(xa, ya, OBJ.');
colorbar;
axis equal tight;
title('OBJECT')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2 -- IDENTIFY EDGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE EDGE MAPS
ECLK = zeros(Nx, Ny);
EOBJ = zeros(Nx, Ny);

% IDENTIFY EDGES
for ny = 2 : Ny - 1
    for nx = 2 : Nx - 1
        if ~CLK(nx,ny)
            A = CLK(nx-1:nx+1, ny-1:ny+1);
            ECLK(nx,ny) = sum(A(:)) > 0;
        end
        
        if OBJ(nx,ny)
            A = OBJ(nx-1:nx+1, ny-1:ny+1);
            EOBJ(nx,ny) = sum(A(:)) < 9;
        end
    end
end

% SHOW EDGES
subplot(1,2,1);
imagesc(xa, ya, (CLK + 2*ECLK).');
colorbar;
axis equal tight;
title('CLOAK + 2*ECLK')

subplot(1,2,2);
imagesc(xa, ya, (OBJ + 2*EOBJ).');
colorbar;
axis equal tight;
title('OBJ + 2*EOBJ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEPS 3 & 4 -- GENERATE BCs FOR SPATIAL TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ORIGINAL COORDINATES
[YA, XA] = meshgrid(ya, xa);

% FORCE MAP
F = ECLK;
XF = F.*XA;
YF = F.*YA;
F = F + EOBJ;


% SHOW MAPS
clf;

subplot(131);
imagesc(xa, ya, F.');
axis equal tight;
colorbar;
title('F - Force Map');

subplot(2,3,2)
imagesc(xa, ya, (F.*XA).');
axis equal tight;
colorbar;
title('XA');

subplot(2,3,3)
imagesc(xa, ya, (F.*YA).');
axis equal tight;
colorbar;
title('YA');

subplot(2,3,5)
imagesc(xa, ya, XF.');
axis equal tight;
colorbar;
title('XF');

subplot(2,3,6)
imagesc(xa, ya, YF.');
axis equal tight;
colorbar;
title('YF');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 5 -- SOLVE LAPLACES EQUATION INSIDE CLOAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IDENTITY MATRIX
M = Nx*Ny;
I = speye(M, M);

% BUILD DERIVATIVE MATRICES
NS = [Nx Ny];
RES = [dx dy];
[DX, D2X, DY, D2Y] = toder(NS, RES);

% INDICES WHERE TO SOLVE LAPLACE EQUATION
PTS = CLK + ECLK - OBJ + EOBJ;
ind = find(PTS(:));

% LAPLACIAN MATRIX
F = diag(sparse(F(:)));
L = D2X + D2Y;
L = F + (I - F)*L;
L = L(ind, ind);

% X TRANSFORM COORDINATES
b       = F*XF(:);
x       = L\b(ind);
XB      = zeros(Nx, Ny);
XB(ind) = x;

% Y TRANSFORM COORDINATES
b       = F*YF(:);
y       = L\b(ind);
YB      = zeros(Nx, Ny);
YB(ind) = y;

% SHOW TRANSFORMED COORDINATES
clf;

subplot(121)
imagesc(xa, ya, XB.');
axis equal tight;
colorbar;
title('XB');

subplot(122)
imagesc(xa, ya, YB.');
axis equal tight;
colorbar;
title('YB');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 6 -- CALCULATE PERMITTIVITY TENSOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE PERMITTIVITY TENSOR
ERxx =  ones(Nx, Ny);   ERxy = zeros(Nx, Ny);   ERxz = zeros(Nx, Ny);
ERyx = zeros(Nx, Ny);   ERyy =  ones(Nx, Ny);   ERyz = zeros(Nx, Ny);
ERzx = zeros(Nx, Ny);   ERzy = zeros(Nx, Ny);   ERzz =  ones(Nx, Ny);

% CALCULATE NUMERICAL DERIVATIVES FOR THE JACOBIAN
DX_XB = reshape(DX*XB(:), Nx, Ny);
DY_XB = reshape(DY*XB(:), Nx, Ny);
DX_YB = reshape(DX*YB(:), Nx, Ny);
DY_YB = reshape(DY*YB(:), Nx, Ny);

% CALCULATE PERMITTIVITY FUNCTION
DEV = CLK - OBJ;
for ny = 1 : Ny
    for nx = 1 : Nx
        if DEV(nx, ny)
            % Get Background Tensor
            ER = [ ERxx(nx, ny) ERxy(nx, ny) ERxz(nx, ny); ...
                   ERyx(nx, ny) ERyy(nx, ny) ERyz(nx, ny); ...
                   ERzx(nx, ny) ERzy(nx, ny) ERzz(nx, ny) ];
               
            % Build Jacobian
            J = [ DX_XB(nx, ny) DY_XB(nx, ny) 0; ...
                  DX_YB(nx, ny) DY_YB(nx, ny) 0; ... 
                  0 0 1 ];
              
            % Transform ER
            J = inv(J);
            ER = J*ER*J.'/det(J);
            
            % Put Values Back into Grid
            ERxx(nx, ny) = ER(1, 1);    ERxy(nx, ny) = ER(1, 2);   ERxz(nx, ny) = ER(1, 3);   
            ERyx(nx, ny) = ER(2, 1);    ERyy(nx, ny) = ER(2, 2);   ERyz(nx, ny) = ER(2, 3);
            ERzx(nx, ny) = ER(3, 1);    ERzy(nx, ny) = ER(3, 2);   ERzz(nx, ny) = ER(3, 3);
        end
    end
end

% VISUALIZE PERMITTIVITY TENSOR
subplot(331);
imagesc(xa, ya, ERxx.');
axis equal tight off;
colorbar;
title('ERxx')

subplot(332);
imagesc(xa, ya, ERxy.');
axis equal tight off;
colorbar;
title('ERxy')

subplot(333);
imagesc(xa, ya, ERxz.');
axis equal tight off;
colorbar;
title('ERxz')

subplot(334);
imagesc(xa, ya, ERyx.');
axis equal tight off;
colorbar;
title('ERyx')

subplot(335);
imagesc(xa, ya, ERyy.');
axis equal tight off;
colorbar;
title('ERyy')

subplot(336);
imagesc(xa, ya, ERyz.');
axis equal tight off;
colorbar;
title('ERyz')


subplot(337);
imagesc(xa, ya, ERzx.');
axis equal tight off;
colorbar;
title('ERzx')

subplot(338);
imagesc(xa, ya, ERzy.');
axis equal tight off;
colorbar;
title('ERzy')

subplot(339);
imagesc(xa, ya, ERzz.');
axis equal tight off;
colorbar;
title('ERzz')




