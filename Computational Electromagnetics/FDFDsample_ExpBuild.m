% demo_gmr
%
% This MATLAB code demonstrates the finite-difference frequency-domain
% method by modeling transmission and reflectance from a guided-mode
% resonance filter.
%
%
% Raymond C. Rumpf
% Prime Research, LC
% Blacksburg, VA
%
% Example code written for short course on
% "Introduction to Optical Simulation Using the Finite-Difference
% Frequency-Domain Method."
% INITIALIZE MATLAB
close all; clc;
clear all;

% UNITS
micrometers = 1;
nanometers = micrometers / 1000;
radians = 1;
degrees = pi/180;

% INITIALIZE FIGURE WINDOW
ss = get(0,'ScreenSize');
fig =figure('Position',[1 0.056*ss(4) ss(3) 0.87*ss(4)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE
LAMBDA = [400:1:700] * nanometers; %wavelength range to simulate
theta = 0 * degrees; %angle of incidence
pol = 'E'; %polarization or mode
% GMR
T = 134 * nanometers; %thickness of GMR grating
L = 314 * nanometers; %period of GMR grating
n1 = 1.00; %refractive index of superstrate
n2 = 1.52; %refractive index of substrate
nL = 2.00; %low refractive index of GMR
nH = 2.10; %high refractive index of GMR
f = 0.5; %duty cycle of GMR grating
% GRID
nmax = max([n1 n2 nL nH]); %maximum refractive index
NRES = 40; %grid resolution parameter
NPML = 20; %size of perfectly matched layer
buf_ylo = 500 * nanometers; %space between GMR and PML
buf_yhi = 500 * nanometers; %space between GMR and PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL GRID RESOLUTION
dx1 = min(LAMBDA) /nmax/NRES; %must resolve smallest wavelength
dy1 = dx1;

dx2 = f*L/NRES; %must resolve finest dimension
dy2 = T/NRES;

dx = min(dx1,dx2); %choose smallest numbers
dy = min(dy1,dy2);

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx); dx = L/nx; %snap x-grid to GMR period
ny = ceil(T/dy); dy = T/ny; %snap y-grid to GMR thickness

% GRID SIZE
Sx = L; %physical size along x
Sy = buf_ylo + T + buf_yhi; %physical size along y
Nx = round(Sx/dx); %grid size along x
Ny = round(Sy/dy) + 2*NPML; %grid size along y
% ENSURE Nx IS ODD FOR FFT
Nx = 2*round(Nx/2) + 1;
dx = Sx/Nx;

% 2X GRID PARAMETERS
Nx2 = 2*Nx; %grid size is twice
Ny2 = 2*Ny;

dx2 = dx/2; %grid spacing is halved
dy2 = dy/2;

% GRID AXES
xa = [0:Nx-1]*dx; xa = xa - mean(xa);
ya = [0:Ny-1]*dy;
xa2 = [0:Nx2-1]*dx2; xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2;

% REPORT GRID SIZE
disp(['Grid size is: Nx = ' num2str(Nx)]);
disp([' Ny = ' num2str(Ny)]);

% CLEAR TEMPORARY VARIABLES
clear dx1 dy1 nx ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD GMR DEVICE ON 2X GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRUCTURE INDICES
nx1 = round((Nx2-f*L/dx2)/2); %position of left of high index
nx2 = round((Nx2+f*L/dx2)/2); %position of right of high index
ny1 = 2*NPML + round(buf_ylo/dy2); %position of top of GMR
ny1 = 2*round(ny1/2) + 1; %make position index odd
ny2 = ny1 + round(T/dy2) - 1; %compute bottom position


% BUILD GMRER2 = er1*ones(Nx2,Ny2);


N2X = n1 * ones(Nx2,Ny2); %initialize entire grid with n1
N2X(:,ny1:ny2) = nL; %fill GMR with nL
N2X(nx1:nx2,ny1:ny2) = nH; %add nH to GMR
N2X(:,ny2+1:Ny2) = n2; %add n2 to superstrate region
% SHOW DEVICE
subplot(141);
imagesc(xa/micrometers,ya/micrometers,N2X');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('DEVICE');
colorbar;
drawnow;

% COMPUTE DIELECTRIC AND MAGNETIC FUNCTIONS
ER2 = N2X.^2; %er = n^2
UR2 = ones(Nx2,Ny2); %no magnetic response
% CLEAR TEMPORARY VARIABLES
clear nx1 nx2 ny1 ny2 N2X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATE GUIDED-MODE RESONANCE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE RECORD VARIABLES
NLAM = length(LAMBDA); %determine how many simulations

REF = zeros(NLAM,1); %initialize reflection record
TRN = zeros(NLAM,1); %initialize transmission record
% ITERATE OVER WAVELENGTH
for nlam = 1 : NLAM

 % Get Next Wavelength
 lam0 = LAMBDA(nlam);

 % Compute Source Vector
 k0 = 2*pi/lam0;
 kinc = k0*n1 * [sin(theta);cos(theta)];

 % Call FDFD2D
 RES2 = [dx2 dy2];
 [R,T,m,F] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,pol);

 % Record Transmission and Reflection Response (%)
 REF(nlam) = 100*sum(R(:)); %add all spatial harmonics
 TRN(nlam) = 100*sum(T(:)); % to compute total power

 % Show Field
 subplot(142);
 imagesc(xa/micrometers,ya/micrometers,real(F)');
 xlabel('x (\mum)');
 ylabel('y (\mum)');
 title('FIELD');
 colorbar;

 % Show Spectra
 subplot(1,4,3:4);
 plot(LAMBDA(1:nlam)/nanometers,REF(1:nlam),'-r'); hold on;
 plot(LAMBDA(1:nlam)/nanometers,TRN(1:nlam),'-b');
 plot(LAMBDA(1:nlam)/nanometers,REF(1:nlam)+TRN(1:nlam),':k'); hold off;
 axis([min(LAMBDA)/nanometers max(LAMBDA)/nanometers 0 105]);
 xlabel('Wavelength (nm)');
 ylabel('% ','Rotation',0);
 title('SPECTRAL RESPONSE');

 drawnow; %update graphics now!
end


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

function [R,T,m,F] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,pol)
% FDFD2D Two-Dimensional Finite-Difference Frequency-Domain
%
% This MATLAB code simulates optical structures using the
% finite-difference frequency-domain method.
%
% INPUT ARGUMENTS
% lam0 is the free space wavelength
% UR2 contains the relative permeability on a 2X grid
% ER2 contains the relative permittivity on a 2X grid
% NPML is the size of the PML on the 1X grid
% RES2 = [dx2 dy2]
% kinc is the indicent wave vector
% pol is the polarization ('E' or 'H')
%
% OUTPUT ARGUMENTS
% R contains diffraction efficiencies of reflected waves
% T contains diffraction efficiencies of transmitted waves
% m contains the indices of the harmonics in R and T
% F is the computed field
%
% Raymond C. Rumpf
% Prime Research, LC
% Blacksburg, VA
%
% Example code written for short course on
% "Introduction to Optical Simulation Using the Finite-Difference
% Frequency-Domain Method."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT AND OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE SIZE OF GRID
[Nx2,Ny2] = size(ER2);
dx2 = RES2(1);
dy2 = RES2(2);

% 1X GRID PARAMETERS
Nx = Nx2/2; dx = 2*dx2;
Ny = Ny2/2; dy = 2*dy2;

% COMPUTE MATRIX SIZE
M = Nx*Ny;

% COMPUTE REFRACTIVE INDEX IN REFLECTION REGION
erref = ER2(:,1); erref = mean(erref(:));
urref = UR2(:,1); urref = mean(urref(:));
nref = sqrt(erref*urref);
if erref<0 && urref<0
 nref = - nref;
end
% COMPUTE REFRACTIVE INDEX IN TRANSMISSION REGION
ertrn = ER2(:,Ny2); ertrn = mean(ertrn(:));
urtrn = UR2(:,Ny2); urtrn = mean(urtrn(:));
ntrn = sqrt(ertrn*urtrn);
if ertrn<0 && urtrn<0
 ntrn = - ntrn;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INCORPORATE PERFECTLY MATCHED LAYER BOUNDARY CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PML PARAMETERS
N0 = 376.73032165; %free space impedance
amax = 3;
cmax = 1;
p = 3;

% INITIALIZE PML TO PROBLEM SPACE
sx = ones(Nx2,Ny2);
sy = ones(Nx2,Ny2);

% COMPUTE FREE SPACE WAVE NUMBERS
k0 = 2*pi/lam0;

% Y PML
N = 2*NPML;
for n = 1 : N
 % compute PML value
 ay = 1 + amax*(n/N)^p;
 cy = cmax*sin(0.5*pi*n/N)^2;
 s = ay*(1-i*cy*N0/k0);
 % incorporate value into PML
 sy(:,N-n+1) = s;
 sy(:,Ny2-N+n) = s;
end
% NGRID2X = [Nx2, Ny2];
% NPML = [0, 0, 20, 20];
% 
% 
% [sx,sy] = calcpml2d(NGRID2X,2*NPML);

% COMPUTE TENSOR COMPONENTS WITH PML
ER2xx = ER2 ./ sx .* sy;
ER2yy = ER2 .* sx ./ sy;
ER2zz = ER2 .* sx .* sy;

UR2xx = UR2 ./ sx .* sy;
UR2yy = UR2 .* sx ./ sy;
UR2zz = UR2 .* sx .* sy;

% OVERLAY MATERIALS ONTO 1X GRID
ERxx = ER2xx(2:2:Nx2,1:2:Ny2);
ERyy = ER2yy(1:2:Nx2,2:2:Ny2);
ERzz = ER2zz(1:2:Nx2,1:2:Ny2);
URxx = UR2xx(1:2:Nx2,2:2:Ny2);
URyy = UR2yy(2:2:Nx2,1:2:Ny2);
URzz = UR2zz(2:2:Nx2,2:2:Ny2);

% CLEAR TEMPORARY VARIABLES
clear N0 amax cmax p sx sy n N ay cy s;
clear UR2 ER2 ER2xx ER2yy ER2zz UR2xx UR2yy UR2zz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FINITE-DIFFERENCE FREQUENCY-DOMAIN ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORM DIAGONAL MATERIAL MATRICES
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% COMPUTE DERIVATIVE OPERATORS
NS = [Nx Ny];
RES = [dx dy];
BC = [-2 -2 0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC,kinc/k0);

% COMPUTE FIELD MATRIX
switch pol
 case 'E',
 A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
 case 'H',
 A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
 otherwise,
 error('Unrecognized polarization.');
end

% COMPUTE SOURCE FIELD
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;
[Y,X] = meshgrid(ya,xa);
fsrc = exp(-i*(kinc(1)*X+kinc(2)*Y));
fsrc = fsrc(:);

% COMPUTE SCATTERED-FIELD MASKING MATRIX
Q = zeros(Nx,Ny);
Q(:,1:NPML+2) = 1;
Q = diag(sparse(Q(:)));

% COMPUTE SOURCE VECTOR
f = (Q*A-A*Q)*fsrc;

% PREPARE MEMORY
clear NS RES BC DEX DEZ DHX DHZ;
clear ya X Y fsrc;
clear ERxx ERyy ERzz URxx URyy URzz;

% COMPUTE FIELD
F = A\f; %backward division is used here!!
F = full(F);
F = reshape(F,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE DIFFRACTION EFFICIENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT REFLECTED AND TRANSMITTED WAVES
Fref = F(:,NPML+1);
Ftrn = F(:,Ny-NPML);

% REMOVE PHASE TILT
p = exp(+i*kinc(1)*xa');
Fref = Fref .* p;
Ftrn = Ftrn .* p;

% COMPUTE SPATIAL HARMONICS
Fref = fftshift(fft(Fref))/Nx;
Ftrn = fftshift(fft(Ftrn))/Nx;

% COMPUTE WAVE VECTOR COMPONENTS OF THE SPATIAL HARMONICS
m = [-floor(Nx/2):floor(Nx/2)]';
kx = kinc(1) - 2*pi*m/(Nx*dx);
kzR = conj( sqrt((k0*nref)^2 - kx.^2) );
kzT = conj( sqrt((k0*ntrn)^2 - kx.^2) );

% COMPUTE DIFFRACTION EFFICIENCY
switch pol
 case 'E',
 R = abs(Fref).^2 .* real(kzR/urref/kinc(2));
 T = abs(Ftrn).^2 .* real(kzT*urref/kinc(2)/urtrn);
 case 'H',
 R = abs(Fref).^2 .* real(kzR/kinc(2));
 T = abs(Ftrn).^2 .* real(kzT*erref/kinc(2)/ertrn);
end
end

% function [sx,sy] = calcpml2d(NGRID,NPML,k0)
% % CALCPML2D Calculate the PML parameters on a 2D grid
% %
% % [sx,sy] = calcpml2d(NGRID,NPML);
% %
% % This MATLAB function calculates the PML parameters sx and sy
% % to absorb outgoing waves on a 2D grid.
% %
% % Input Arguments
% % =================
% % NGRID Array containing the number of points in the grid
% % = [ Nx Ny ]
% % NPML Array containing the size of the PML at each boundary
% % = [ Nxlo Nxhi Nylo Nyhi ]
% %
% % Output Arguments
% % =================
% % sx,sy 2D arrays containing the PML parameters on a 2D grid
%     
%     a_max = 3;
%     p = 3;
%     sigmaprime_max = 1;
%     eta0 = 376.73032165;
% 
%     sx = ones(NGRID);
%     sy = sx ;
%     
%     for nx = 1:NPML(1)
%         sx(NPML(1)-nx+1,:) = (1 + a_max*(nx/NPML(1))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*nx)/(2*NPML(1))))^2));
% 
%     end
% 
%     for nx = 1 : NPML(2)
%         sx(NGRID(1)-(2)+nx,:) = (1 + a_max*(nx/NPML(2))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*nx)/(2*NPML(2))))^2));
% 
%     end
% 
%     for ny = 1:NPML(3)
%         sy(:,NPML(3)-ny+1) = (1 + a_max*(ny/NPML(3))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*ny)/(2*NPML(3))))^2));
% 
%     end
% 
%     for ny = 1 : NPML(4)
%         sy(:,NGRID(2)-NPML(4)+ny) = (1 + a_max*(ny/NPML(4))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*ny)/(2*NPML(4))))^2));
% 
%     end
% end
