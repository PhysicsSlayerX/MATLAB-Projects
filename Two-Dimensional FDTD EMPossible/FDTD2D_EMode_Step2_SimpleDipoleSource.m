% FDTD2D_EMode_Step2_SimpleDipoleSource.m
% Inject a simple dipole source

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

% CONSTANTS
e0 = 8.85418782e-12 * 1/meters;
u0 = 1.25663706e-6 * 1/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
f0      = 1.0 * gigahertz;
lam0    = c0/f0;

% GRID PARAMETERS
NRES    = 20;
Sx      = 10*lam0;
Sy      = 10*lam0;
nmax    = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% COMPUTE GRID SIZE
Nx = ceil(Sx/dx);
Sx = Nx*dx;
Ny = ceil(Sy/dy);
Sy = Ny*dy;

% GRID AXES
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO VACUUM

ERzz = ones(Nx,Ny);
URxx = ones(Nx,Ny);
URyy = ones(Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE STABLE TIME STEP
dmin    = min([dx dy]);
dt      = dmin/(2*c0);

% CALCULATE STABLE TIME STEP
nx_src  = round(0.3*Nx);
ny_src  = round(0.4*Ny);
tau     = 0.5/f0;
t0      = 3*tau;

% CALCULATE NUMBER OF TIME STEPS
tprop   = nmax*Sy/c0;
t       = 2*t0 + 2*tprop;
STEPS   = ceil(t/dt);

% CALCULATE GAUSSIAN PULSE SOURCE
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

% INITIALIZE FIELDS TO ZERO
Bx = zeros(Nx, Ny);
By = zeros(Nx, Ny);
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);

% INITIALIZE DERIVATIVE ARRAYS
dEzx = zeros(Nx, Ny);
dEzy = zeros(Nx, Ny);
dHxy = zeros(Nx, Ny);
dHyx = zeros(Nx, Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN FDTD LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATE OVER TIME
%
for T = 1 : STEPS
   
    % CALCULATE DERIVATIVES OF E
        % dEzy
        for nx = 1: Nx
            for ny  = 1 : Ny-1
                dEzy(nx,ny) = ( Ez(nx, ny+1) - Ez(nx,ny) )/dy;
            end
            dEzy(nx,Ny) = ( 0 - Ez(nx,Ny) )/dy;

        end
        % dEzx
        for ny = 1 : Ny
            for nx = 1 : Nx-1
                dEzx(nx,ny) = ( Ez(nx+1,ny) - Ez(nx,ny) )/dx;
            end
            dEzx(Nx,ny) = ( 0 - Ez(Nx,ny) )/dx;
        end
      
    % UPDATE B FROM E
    Bx = Bx - M*(dEzy);
    By = By - M*(-dEzx);
    
    % UPDATE H FROM B
    Hx = Bx./URxx;
    Hy = By./URyy;
        
    % CALCULATE DERIVATIVES OF H
        % dHyx
        for ny = 1 : Ny
            dHyx(1,ny) = ( Hy(1,ny) - 0 )/dx;
            for nx = 2 : Nx
                dHyx(nx,ny) = ( Hy(nx,ny) - Hy(nx-1,ny) )/dx;
            end
        end
        % dHxy
        for nx = 1 : Nx
            dHxy(nx,1) = ( Hx(nx,1) - 0 )/dy;
            for ny = 2 : Ny
                dHxy(nx,ny) = ( Hx(nx,ny) - Hx(nx,ny-1) )/dy;
            end
        end
        
    % UPDATE D FROM H
    Dz = Dz + M*(dHyx - dHxy);
    
    % UPDATE E FROM D
    Ez = Dz./ERzz;
    
    % INJECT SIMPLE DIPOLE SOURCE
    Ez(nx_src,ny_src) = Ez(nx_src, ny_src) + gsrc(T);
    
    % SHOW FIELDS
    if mod(T, 20) == 0
        subplot(1,3,1);
        imagesc(xa,ya,Ez.');
        axis equal tight;
        colorbar;
        title('Ez');
        caxis(0.002*[-1 +1]);
        
        subplot(1,3,2);
        imagesc(xa,ya,Hx.');
        axis equal tight;
        colorbar;
        title('Hx');
        caxis(0.002*[-1 +1]);
        
        subplot(1,3,3);
        imagesc(xa,ya,Hy.');
        axis equal tight;
        colorbar;
        title('Hy');
        caxis(0.002*[-1 +1]);
        
        drawnow;
    end
    
end









