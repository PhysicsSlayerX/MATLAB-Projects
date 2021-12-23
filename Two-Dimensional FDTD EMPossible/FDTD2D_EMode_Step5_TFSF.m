% FDTD2D_EMode_Step6_TFSF.m
% Incorporate Total Field/Scattered Field

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
pml_kx   = 1;
pml_ky   = 1;
pml_ax   = 1e-10;
pml_ay   = 1e-10;
pml_Npml = 3;
pml_R0   = 1e-8;

% GRID PARAMETERS
NRES    = 20;
Sx      = 10*lam0;
Sy      = 10*lam0;
NPML    = [20 21 22 23];
nmax    = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% COMPUTE GRID SIZE
Nx = NPML(1) + ceil(Sx/dx) + NPML(2);
Sx = Nx*dx;
Ny = NPML(3) + ceil(Sy/dy) + NPML(4);
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

% CALCULATE GAUSSIAN PULSE PARAMETERS
nx_src  = NPML(1) + 1;
ny_src  = NPML(3) + 1;
ETAsrc  = sqrt(URxx(nx_src,ny_src)/ERzz(nx_src, ny_src)); 
nsrc  = sqrt(URxx(nx_src,ny_src)*ERzz(nx_src, ny_src)); 
tau     = 0.5/f0;
t0      = nsrc*sqrt((Sx/2)^2 + (Sy/2)^2)/c0;

% CALCULATE NUMBER OF TIME STEPS
tprop   = nmax*Sy/c0;
t       = 2*t0 + 2*tprop;
STEPS   = ceil(t/dt);

% DIRECTION COSINES
gx = cos(theta);
gy = sin(theta);

% COMPENSATE FOR DISPERSION
k0 = 2*pi*f0/c0;
f  = c0*dt/sin(pi*f0*dt)*sqrt((sin(k0*gx*dx/2)/dx)^2 ...
                             +(sin(k0*gy*dy/2)/dy)^2);

ERzz = f*ERzz;
URxx = f*URxx;
URyy = f*URyy;
       
% Q
[Y, X] = meshgrid(ya - mean(ya), xa - mean(xa));
r = 0.45*(Sx - NPML(1)*dx - NPML(2)*dx);
Q = (X.^2 + Y.^2) >= r^2;

% INITIALIZE TF/SF INTERFACE INDEX ARRAYS
is2t = zeros(Nx, Ny);
it2s = zeros(Nx, Ny);
js2t = zeros(Nx, Ny);
jt2s = zeros(Nx, Ny);

isft = zeros(Nx, Ny);
itfs = zeros(Nx, Ny);
jsft = zeros(Nx, Ny);
jtfs = zeros(Nx, Ny);

% DETECT TF/SF INTERFACE CELLS
is2t(1:Nx-1, :)     = Q(1:Nx-1, :) > Q(2:Nx,:);
it2s(1:Nx-1, :)     = Q(1:Nx-1, :) < Q(2:Nx,:);
js2t(:,1:Ny-1)      = Q(:,1:Ny-1) > Q(:, 2:Ny);
jt2s(:,1:Ny-1)      = Q(:,1:Ny-1) < Q(:, 2:Ny);

isft(2:Nx,:)        = Q(2:Nx,:) > Q(1:Nx-1,:);
itfs(2:Nx,:)        = Q(2:Nx,:) < Q(1:Nx-1,:);
jsft(:,2:Ny)        = Q(:,2:Ny) > Q(:,1:Ny-1);
jtfs(:,2:Ny)        = Q(:,2:Ny) < Q(:,1:Ny-1);

% TABULATE LINEAR ARRAY INDICES
is2t = find(is2t(:));
it2s = find(it2s(:));
js2t = find(js2t(:));
jt2s = find(jt2s(:));

isft = find(isft(:));
itfs = find(itfs(:));
jsft = find(jsft(:));
jtfs = find(jtfs(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PML CONDUCTIVITY ON 2X GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE PML CONDUCTIVITIES
sigx2 = zeros(Nx2,Ny2);
sigy2 = zeros(Nx2,Ny2);

% ADD XLO CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(1)*dx2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx2(nx1,:) = sigmax*(nx/2/NPML(1))^pml_Npml;
end

% ADD XHI CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(2)*dx2);
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx2(nx1,:) = sigmax*(nx/2/NPML(2))^pml_Npml;
end

% ADD YLO CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(3)*dy2);
for ny = 1 : 2*NPML(3)
    ny1 = 2*NPML(3) - ny + 1;
    sigy2(:,ny1) = sigmax*(ny/2/NPML(3))^pml_Npml;
end

% ADD YHI CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(4)*dy2);
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy2(:,ny1) = sigmax*(ny/2/NPML(4))^pml_Npml;
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M
M = c0*dt;

% CALCULATE UPDATE COEFFICIENTS FOR Bx
sigx = sigx2(1:2:Nx2, 2:2:Ny2);
sigy = sigy2(1:2:Nx2, 2:2:Ny2);
bBxy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cBxy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bBxy - 1);

% CALCULATE UPDATE COEFFICIENTS FOR By
sigx = sigx2(2:2:Nx2, 1:2:Ny2);
sigy = sigy2(2:2:Nx2, 1:2:Ny2);
bByx = exp(-(sigx/pml_kx + pml_ax)*dt/e0);
cByx = sigx./(sigx*pml_kx + pml_ax*pml_kx^2).*(bByx - 1);

% CALCULATE UPDATE COEFFICIENTS FOR Dz
sigx = sigx2(1:2:Nx2,1:2:Ny2);
sigy = sigy2(1:2:Nx2,1:2:Ny2);
bDzx = exp(-(sigx/pml_kx + pml_ax)*dt/e0);
cDzx = sigx./(sigx*pml_kx + pml_ax*pml_kx^2).*(bDzx - 1);
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
psiBx_ylo = zeros(Nx,NPML(3));
psiBx_yhi = zeros(Nx,NPML(4));

psiBy_xlo = zeros(NPML(1), Ny);
psiBy_xhi = zeros(NPML(2), Ny);

psiDz_xlo = zeros(NPML(1), Ny);
psiDz_xhi = zeros(NPML(2), Ny);
psiDz_ylo = zeros(Nx,NPML(3));
psiDz_yhi = zeros(Nx,NPML(4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN FDTD LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATE OVER TIME
%
for T = 1 : STEPS
   
    % CALCULATE SOURCE FIELDS
    Ezsrc =            exp(-(((T + 0  )*dt - t0 - (gx*(X       ) + gy*(Y       ))*nsrc/c0)/tau).^2);
    Hxsrc = +gy/ETAsrc*exp(-(((T + 0.5)*dt - t0 - (gx*(X       ) + gy*(Y + dy/2))*nsrc/c0)/tau).^2);
    Hysrc = -gx/ETAsrc*exp(-(((T + 0.5)*dt - t0 - (gx*(X + dx/2) + gy*(Y       ))*nsrc/c0)/tau).^2);
    
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
    dEzx(is2t) = dEzx(is2t) - Ezsrc(is2t+1)/dx;
    dEzx(it2s) = dEzx(it2s) + Ezsrc(it2s+1)/dx;
    dEzy(js2t) = dEzy(js2t) - Ezsrc(js2t+Nx)/dy;
    dEzy(jt2s) = dEzy(jt2s) + Ezsrc(jt2s+Nx)/dy;
      
    % UPDATE CONVOLUTIONS FOR B FIELD UPDATES
    psiBx_ylo = bBxy(:,1:NPML(3)).*psiBx_ylo ...
              + cBxy(:,1:NPML(3)).*dEzy(:,1:NPML(3));
    psiBx_yhi = bBxy(:,Ny-NPML(4)+1:Ny).*psiBx_yhi ...
              + cBxy(:,Ny-NPML(4)+1:Ny).*dEzy(:,Ny-NPML(4)+1:Ny);        
        
    psiBy_xlo = bByx(1:NPML(1),:).*psiBy_xlo ...
              + cByx(1:NPML(1),:).*dEzx(1:NPML(1),:);    
    psiBy_xhi = bByx(Nx-NPML(2)+1:Nx,:).*psiBy_xhi ...
              + cByx(Nx-NPML(2)+1:Nx,:).*dEzx(Nx-NPML(2)+1:Nx,:);    
            
    % UPDATE B FROM E
    Bx                      = Bx - M*(dEzy/pml_ky);
    Bx(:,1:NPML(3))         = Bx(:,1:NPML(3)) - M*psiBx_ylo;
    Bx(:,Ny-NPML(4)+1:Ny)   = Bx(:,Ny-NPML(4)+1:Ny) - M*psiBx_yhi;
    
    By                      = By - M*(-dEzx/pml_kx);
    By(1:NPML(1),:)         = By(1:NPML(1),:) + M*psiBy_xlo;
    By(Nx-NPML(2)+1:Nx,:)   = By(Nx-NPML(2)+1:Nx,:) + M*psiBy_xhi;
    
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
    dHyx(itfs) = dHyx(itfs) - Hysrc(itfs-1)/dx;
    dHyx(isft) = dHyx(isft) + Hysrc(isft-1)/dx;    
    dHxy(jtfs) = dHxy(jtfs) - Hxsrc(jtfs-Nx)/dy;
    dHxy(jsft) = dHxy(jsft) + Hxsrc(jsft-Nx)/dy;
    
    % UPDATE CONVOLUITIONS FOR D FIELD UPDATE
    psiDz_xlo = bDzx(1:NPML(1),:).*psiDz_xlo ...
              + cDzx(1:NPML(1),:).*dHyx(1:NPML(1),:);
    psiDz_xhi = bDzx(Nx - NPML(2)+1:Nx,:).*psiDz_xhi ...
              + cDzx(Nx - NPML(2)+1:Nx,:).*dHyx(Nx - NPML(2)+1:Nx,:);
    psiDz_ylo = bDzy(:,1:NPML(3)).*psiDz_ylo ...
              + cDzy(:,1:NPML(3)).*dHxy(:,1:NPML(3));      
    psiDz_yhi = bDzy(:,Ny-NPML(4)+1:Ny).*psiDz_yhi ...
              + cDzy(:,Ny-NPML(4)+1:Ny).*dHxy(:,Ny-NPML(4)+1:Ny);      
    
          
    % UPDATE D FROM H
    Dz                      = Dz + M*(dHyx/pml_kx - dHxy/pml_ky);
    Dz(1:NPML(1),:)         = Dz(1:NPML(1),:) + M*psiDz_xlo;
    Dz(Nx - NPML(2)+1:Nx,:) = Dz(Nx - NPML(2)+1:Nx,:) + M*psiDz_xhi;
    Dz(:,1:NPML(3))         = Dz(:,1:NPML(3)) - M*psiDz_ylo;
    Dz(:,Ny-NPML(4)+1:Ny)   = Dz(:,Ny-NPML(4)+1:Ny) - M*psiDz_yhi;
    
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









