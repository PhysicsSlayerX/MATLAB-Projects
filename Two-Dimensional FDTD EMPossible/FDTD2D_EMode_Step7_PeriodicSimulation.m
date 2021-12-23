% FDTD2D_EMode_Step7_PeriodicSimulation.m
% Simplify Scattering Code for Periodic Structures

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
theta   = 0*degrees;

% PML PARAMETERS
pml_ky   = 1;
pml_ay   = 1e-10;
pml_Npml = 3;
pml_R0   = 1e-8;

% GRID PARAMETERS
NRES    = 20;
Sx      = 3*lam0;
SPACER  = lam0 * [2 2];
NPML    = [20 20];
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

Sy = SPACER(1) + 0 + SPACER(2);
Ny = NPML(1) + ceil(Sy/dy) + NPML(2);
Sy = Ny*dy;

% 2X GRID PARAMETERS
Nx2 = 2*Nx;     dx2 = dx/2;
Ny2 = 2*Ny;     dy2 = dy/2;

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

% POSITION OF SOURCE
ny_src  = round(0.5*Ny);

% GAUSSIAN PULSE PARAMETERS
tau     = 0.5/f0;
t0      = 3*tau;

% CALCULATE NUMBER OF TIME STEPS
tprop   = nmax*Sy/c0;
t       = 2*t0 + 2*tprop;
STEPS   = ceil(t/dt);

% COMPUTE GAUSSIAN PULSE SOURCE
ETAsrc  = sqrt(URxx(1,ny_src)/ERzz(1, ny_src)); 
nsrc    = sqrt(URxx(1,ny_src)*ERzz(1, ny_src));
delt    = 0.5*dt +  0.5*nsrc*dy/c0;
t       = [0:STEPS-1]*dt;
Ezsrc   = exp(-((t - t0)/tau).^2);
Hxsrc   = (1/ETAsrc)*exp(-((t - t0 + delt)/tau).^2);

% COMPENSATE FOR DISPERSION
k0   = 2*pi*f0/c0;
f    = c0*dt/sin(pi*f0*dt)*sin(k0*dy/2)/dy;
ERzz = f*ERzz;
URxx = f*URxx;
URyy = f*URyy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PML CONDUCTIVITY ON 2X GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE PML CONDUCTIVITIES
sigy2 = zeros(Nx2,Ny2);


% ADD YLO CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(1)*dy2);
for ny = 1 : 2*NPML(1)
    ny1 = 2*NPML(1) - ny + 1;
    sigy2(:,ny1) = sigmax*(ny/2/NPML(1))^pml_Npml;
end

% ADD YHI CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(2)*dy2);
for ny = 1 : 2*NPML(2)
    ny1 = Ny2 - 2*NPML(2) + ny;
    sigy2(:,ny1) = sigmax*(ny/2/NPML(2))^pml_Npml;
end


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M
M = c0*dt;

% CALCULATE UPDATE COEFFICIENTS FOR Bx
sigy = sigy2(1:2:Nx2, 2:2:Ny2);
bBxy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cBxy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bBxy - 1);

% CALCULATE UPDATE COEFFICIENTS FOR Dz
sigy = sigy2(1:2:Nx2,1:2:Ny2);
bDzy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cDzy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bDzy - 1);

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

% INITIALIZE CONVOLUTIONS
psiBx_ylo = zeros(Nx,NPML(1));
psiBx_yhi = zeros(Nx,NPML(2));
psiDz_ylo = zeros(Nx,NPML(1));
psiDz_yhi = zeros(Nx,NPML(2));

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
            dEzy(nx,Ny) = ( Ez(nx,1) - Ez(nx,Ny) )/dy;

        end
        % dEzx
        for ny = 1 : Ny
            for nx = 1 : Nx-1
                dEzx(nx,ny) = ( Ez(nx+1,ny) - Ez(nx,ny) )/dx;
            end
            dEzx(Nx,ny) = ( Ez(1,ny) - Ez(Nx,ny) )/dx;
        end
    % INCORPORATE TF/SF CORRRECTION
    dEzy(:, ny_src-1) = dEzy(:, ny_src-1) - Ezsrc(T)/dy;
      
    % UPDATE CONVOLUTIONS FOR B FIELD UPDATES
    psiBx_ylo = bBxy(:,1:NPML(1)).*psiBx_ylo ...
              + cBxy(:,1:NPML(1)).*dEzy(:,1:NPML(2));
    psiBx_yhi = bBxy(:,Ny-NPML(2)+1:Ny).*psiBx_yhi ...
              + cBxy(:,Ny-NPML(2)+1:Ny).*dEzy(:,Ny-NPML(2)+1:Ny);         
            
    % UPDATE B FROM E
    Bx                      = Bx - M*(dEzy/pml_ky);
    Bx(:,1:NPML(1))         = Bx(:,1:NPML(1)) - M*psiBx_ylo;
    Bx(:,Ny-NPML(2)+1:Ny)   = Bx(:,Ny-NPML(2)+1:Ny) - M*psiBx_yhi;
    
    By                      = By - M*(-dEzx);
    
    % UPDATE H FROM B
    Hx = Bx./URxx;
    Hy = By./URyy;
        
    % CALCULATE DERIVATIVES OF H
        % dHyx
        for ny = 1 : Ny
            dHyx(1,ny) = ( Hy(1,ny) - Hy(Nx,ny) )/dx;
            for nx = 2 : Nx
                dHyx(nx,ny) = ( Hy(nx,ny) - Hy(nx-1,ny) )/dx;
            end
        end
        % dHxy
        for nx = 1 : Nx
            dHxy(nx,1) = ( Hx(nx,1) - Hx(nx,Ny) )/dy;
            for ny = 2 : Ny
                dHxy(nx,ny) = ( Hx(nx,ny) - Hx(nx,ny-1) )/dy;
            end
        end
        
    % INCORPORATE TF/SF CORRECTIONS
    dHxy(:,ny_src) = dHxy(:,ny_src) - Hxsrc(T)/dy;
    
    
    % UPDATE CONVOLUITIONS FOR D FIELD UPDATE
    psiDz_ylo = bDzy(:,1:NPML(1)).*psiDz_ylo ...
              + cDzy(:,1:NPML(1)).*dHxy(:,1:NPML(2));      
    psiDz_yhi = bDzy(:,Ny-NPML(2)+1:Ny).*psiDz_yhi ...
              + cDzy(:,Ny-NPML(2)+1:Ny).*dHxy(:,Ny-NPML(2)+1:Ny);      
          
    % UPDATE D FROM H
    Dz                      = Dz + M*(dHyx - dHxy/pml_ky);
    Dz(:,1:NPML(1))         = Dz(:,1:NPML(1)) - M*psiDz_ylo;
    Dz(:,Ny-NPML(2)+1:Ny)   = Dz(:,Ny-NPML(2)+1:Ny) - M*psiDz_yhi;
    
    % UPDATE E FROM D
    Ez = Dz./ERzz;
    
    % SHOW FIELDS
    if mod(T, 20) == 0
        subplot(1,3,1);
        imagesc(xa,ya,Ez.');
        axis equal tight;
        colorbar;
        title('Ez');
        caxis([-1 +1]);
        
        subplot(1,3,2);
        imagesc(xa,ya,Hx.');
        axis equal tight;
        colorbar;
        title('Hx');
        caxis([-1 +1]);
        
        subplot(1,3,3);
        imagesc(xa,ya,Hy.');
        axis equal tight;
        colorbar;
        title('Hy');
        caxis([-1 +1]);
        
        drawnow;
    end
    
end









