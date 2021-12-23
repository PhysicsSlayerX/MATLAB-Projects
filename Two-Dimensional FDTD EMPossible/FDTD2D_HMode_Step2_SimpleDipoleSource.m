% FDTD2D_HMode_Step2_SimpleDipoleSource.m
% Inject a Simple Dipole Source

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees         = pi/180;
meters          = 1;
centimeters     = 1e-2 * meters;
millimmeters    = 1e-3 * meters;
micrometers     = 1e-6 * meters;
nanometers      = 1e-9 * meters;
inches          = 2.54 * centimeters;
feet            = 12 * inches;
seconds         = 1;
hertz           = 1/seconds;
kilohertz       = 1e3 * hertz;
megahertz       = 1e6 * hertz;
gigahertz       = 1e9 * hertz;
terahertz       = 1e12 * hertz;
petahertz       = 1e15 * hertz;

% CONSTANTS
e0 = 8.85418782e-12 * 1/meters;
u0 = 1.25663706e-6 * 1/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
f0 = 1.0 * gigahertz;
lam0 = c0/f0;

% GRID PARAMETERS
NRES     = 20;
Sx       = 10 * lam0;
Sy       = 10 * lam0;
nmax     = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% COMPUTER GRID SIZE
Nx = ceil(Sx/dx);
Sx = Nx*dx;
Ny = ceil(Sy/dy);
Sy = Ny*dy;

% GRID AXES
xa = [0 : Nx - 1]*dx;
ya = [0 : Ny - 1]*dy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO VACUUM
ERxx = ones(Nx,Ny);
ERyy = ones(Nx,Ny);
URzz = ones(Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE STABLE TIME STEP
dmin     = min([dx dy]);
dt       = dmin/(2*c0);

% CALCULATE GAUSSIAN PULSE PARAMETERS
nx_src   = round(0.3*Nx);
ny_src   = round(0.4*Ny);
tau      = 0.5/f0;
t0       = 3*tau;

% CALCULATE NUMBER OF TIME STEPS
tprop = nmax*Sy/c0;
t     = 2*t0 + 2*tprop;
STEPS = ceil(t/dt);

% CALCULATE GAUSSIAN SOURCE
t = [0:STEPS-1]*dt;
gsrc = exp(-((t - t0)/tau).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M
M = c0*dt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE FDTD TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE FIELDS
Bz = zeros(Nx, Ny);
Hz = zeros(Nx, Ny);
Dx = zeros(Nx, Ny);
Dy = zeros(Nx, Ny);
Ex = zeros(Nx, Ny);
Ey = zeros(Nx, Ny);

% INITIALIZE DERIVATIVE ARRAYS
dHzx = zeros(Nx, Ny);
dHzy = zeros(Nx, Ny);
dExy = zeros(Nx, Ny);
dEyx = zeros(Nx, Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN FDTD LOOP FOR H MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATE OVER TIME
%

for T = 1 : STEPS
    
    % CALCULATE DERIVATIVES OF E
        % dEyx
        for ny = 1 : Ny
            for nx = 1 : Nx-1
                dEyx(nx, ny) = ( Ey(nx+1, ny) - Ey(nx,ny) )/dx;
            end
            dEyx(Nx, ny) = ( 0 - Ey(nx,ny) )/dx;
        end
        
        % dExy
        for nx = 1 : Nx
            for ny = 1 : Ny - 1
                dExy(nx, ny) = ( Ex(nx, ny+1) - Ex(nx,ny) ) /dy; 
            end
            dExy(nx, Ny) = ( 0 - Ex(nx,Ny) ) /dy; 
        end
        
     % UPDATE B FROM E
     Bz = Bz - M*(dEyx - dExy);
     
     % UPDATE H FROM B
     Hz = Bz./URzz;
     
     % INJECT SIMPLE DIPOLE SOURCE
     Hz(nx_src, ny_src) = Hz(nx_src, ny_src) + gsrc(T);
     
     % CALCULATE DERIVATIVES OF H
        % dHzy
        for nx = 1 : Nx
            dHzy(nx, 1) = ( Hz(nx,1) - 0 )/dy;
            for ny = 2 : Ny
                dHzy(nx, ny) = ( Hz(nx,ny) - Hz(nx, ny-1) )/dy;
            end
        end
        
        % dHzx
        for ny = 1 : Ny
            dHzx(1, ny) = ( Hz(1,ny) - 0 ) / dx;
            for nx = 2 : Nx
                dHzx(nx, ny) = ( Hz(nx,ny) - Hz(nx-1, ny) )/ dx;
            end
        end
        
     % UPDATE D FROM H
     Dx = Dx + M*(dHzy);
     Dy = Dy + M*(-dHzx);
     
     % UPDATE E FROM D
     Ex = Dx./ERxx;
     Ey = Dy./ERyy;
     
     % SHOW FIELDS
     if mod(T, 20) == 0
         
         subplot(1,3,1);
         imagesc(xa, ya, Hz.')
         axis equal tight
         colorbar;
         title('Hz')
         caxis(0.002*[-1 +1]);
         
         subplot(1,3,2);
         imagesc(xa, ya, Ex.')
         axis equal tight
         colorbar;
         title('Ex')
         caxis(0.002*[-1 +1]);
         
         subplot(1,3,3);
         imagesc(xa, ya, Ey.')
         axis equal tight
         colorbar;
         title('Ey')
         caxis(0.002*[-1 +1]);
         
         drawnow;
     end
end
















