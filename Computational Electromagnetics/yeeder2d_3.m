

degrees = pi/180;

k0 = 2*pi;
theta = 20*degrees;
kinc = k0*[sin(theta);cos(theta)];


NGRID = [12 10];
RES = [1 1];
BC = [-2 -2 -2 -2];


[DEX,DEY,DHX,DHY] = yeeder2d(NGRID,RES,BC,kinc);

dex = full(DEX);
dey = full(DEY);
dhx = full(DHX);
dhy = full(DHY);

xa = [-Nx/2:Nx/2]*-dx;
ya = [0:Ny-1]*-dy;

subplot (1,4,1);
g = imagesc(xa,ya, real(dex).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DX');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));
subplot (1,4,2);
g = imagesc(xa,ya, real(dey).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DY');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));

subplot (1,4,3);
g = imagesc(xa,ya, real(dhx).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DX');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));
subplot (1,4,4);
g = imagesc(xa,ya, real(dhy).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DY');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));

set(gcf, 'Position',  [300, 200, 1500, 700])

function [DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC,kinc);
% YEEDER2D Yee Grid Derivative Operators on a 2D Grid
%
% [DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC,kinc);
%
% Input Arguments
% =================
% NS [Nx Ny] 1X grid size
% RES [dx dy] 1X grid resolution
% BC [xlo xhi ylo yhi] boundary conditions
% -2: pseudo-periodic (requires kinc)
% -1: periodic
% 0: Dirichlet
% kinc [kx ky] incident wave vector
% This argument is only needed for pseudo-periodic boundaries.
%
% Note: For normalized grids, use dx=k0*dx and kinc=kinc/k0
%
% Output Arguments
% =================
%
% Ey(i+1,j) - Ey(i,j) Ex(i,j+1) - Ex(i,j)
% DEX*Ex = ------------------- DEY*Ey = -------------------
% dx dy
%
% Hy(i,j) - Hy(i-1,j) Hx(i,j) - Hx(i,j-1)
% DHX*Hx = ------------------- DHY*Hy = -------------------
% dx dy
%
% Raymond C. Rumpf
% Prime Research, LC
% Blacksburg, VA
%
% Example code written for short course on
% "Introduction to Optical Simulation Using the Finite-Difference
% Frequency-Domain Method."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFY INPUT/OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFY NUMBER OF INPUT ARGUMENTS
error(nargchk(3,4,nargin));

% VERIFY NUMBER OF OUTPUT ARGUMENTS
error(nargchk(1,4,nargout));

% EXTRACT GRID PARAMETERS
Nx = NS(1); dx = RES(1);
Ny = NS(2); dy = RES(2);

% DETERMINE MATRIX SIZE
M = Nx*Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MATRIX
DEX = sparse(M,M);

% PLACE MAIN DIAGONALS
DEX = spdiags(-ones(M,1),0,DEX);
DEX = spdiags(+ones(M,1),+1,DEX);

% CORRECT BOUNDARY TERMS (DEFAULT TO DIRICHLET)
for ny = 1 : Ny-1
 neq = Nx*(ny-1) + Nx;
 DEX(neq,neq+1) = 0;
end

% HANDLE BOUNDARY CONDITIONS ON XHI SIDE
switch BC(2)
 case -2,
 dpx = exp(-i*kinc(1)*Nx*dx);
 for ny = 1 : Ny
 neq = Nx*(ny-1) + Nx;
 nv = Nx*(ny-1) + 1;
 DEX(neq,nv) = +dpx;
 end
 case -1,
 for ny = 1 : Ny
 neq = Nx*(ny-1) + Nx;
 nv = Nx*(ny-1) + 1;
 DEX(neq,nv) = +1;
 end
 case 0, %Dirichlet
 otherwise,
 error('Unrecognized x-high boundary condition.');
end
% FINISH COMPUTATION
DEX = DEX / dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MATRIX
DEY = sparse(M,M);

% PLACE MAIN DIAGONALS
DEY = spdiags(-ones(M,1),0,DEY);
DEY = spdiags(+ones(M,1),+Nx,DEY);

% HANDLE BOUNDARY CONDITIONS ON YHI SIDE
switch BC(4)
 case -2,
 dpy = exp(-i*kinc(2)*Ny*dy);
 for nx = 1 : Nx
 neq = Nx*(Ny-1) + nx;
 nv = nx;
 DEY(neq,nv) = +dpy;
 end
 case -1,
 for nx = 1 : Nx
 neq = Nx*(Ny-1) + nx;
 nv = nx;
 DEY(neq,nv) = +1;
 end
 case 0,
 otherwise,
 error('Unrecognized y-high boundary condition.');
end
% FINISH COMPUTATION
DEY = DEY / dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DHX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MATRIX
DHX = sparse(M,M);

% PLACE MAIN DIAGONALS
DHX = spdiags(+ones(M,1),0,DHX);
DHX = spdiags(-ones(M,1),-1,DHX);

% CORRECT BOUNDARY TERMS (DEFAULT TO DIRICHLET)
for ny = 2 : Ny
 neq = Nx*(ny-1) + 1;

 DHX(neq,neq-1) = 0;
end
% HANDLE BOUNDARY CONDITIONS ON XLOW SIDE
switch BC(1)
 case -2,
 dpx = exp(+i*kinc(1)*Nx*dx);
 for ny = 1 : Ny
 neq = Nx*(ny-1) + 1;
 nv = Nx*(ny-1) + Nx;
 DHX(neq,nv) = -dpx;
 end
 case -1,
 for ny = 1 : Ny
 neq = Nx*(ny-1) + 1;
 nv = Nx*(ny-1) + Nx;
 DHX(neq,nv) = -1;
 end
 case 0,
 otherwise,
 error('Unrecognized x-low boundary condition.');
end
% FINISH COMPUTATION
DHX = DHX / dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MATRIX
DHY = sparse(M,M);

% PLACE MAIN DIAGONALS
DHY = spdiags(+ones(M,1),0,DHY);
DHY = spdiags(-ones(M,1),-Nx,DHY);

% HANDLE BOUNDARY CONDITIONS ON YLOW SIDE
switch BC(3)
 case -2,
 dpy = exp(+i*kinc(2)*Ny*dy);
 for nx = 1 : Nx
 neq = nx;
 nv = Nx*(Ny-1) + nx;
 DHY(neq,nv) = -dpy;
 end
 case -1,
 for nx = 1 : Nx
 neq = nx;
 nv = Nx*(Ny-1) + nx;
 DHY(neq,nv) = -1;
 end
 case 0,
 otherwise,
 error('Unrecognized y-low boundary condition.');
end
% FINISH COMPUTATION
DHY = DHY / dy;
end