% FDTD2D_HMode_Step5_TFSF.m
% Incorporate Total Field/Scattered Field Source
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
f0      = 1.0 * gigahertz;
lam0    = c0/f0;
theta   = 30*degrees;

% PML PARAMETERS
pml_kx      = 1;
pml_ky      = 1;
pml_ax      = 1e-10;
pml_ay      = 1e-10;
pml_Npml    = 3;
pml_R0      = 1e-8;

% GRID PARAMETERS
NRES     = 20;
Sx       = 10 * lam0;
Sy       = 10 * lam0;
NPML     = [ 20 21 22 23 ];
nmax     = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% COMPUTER GRID SIZE
Nx = NPML(1) + ceil(Sx/dx) + NPML(2);
Sx = Nx*dx;
Ny = NPML(3) + ceil(Sy/dy) + NPML(4);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2 * Nx; dx2 = dx/2;
Ny2 = 2 * Ny; dy2 = dy/2;

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
nx_src   = NPML(1) + 1;
ny_src   = NPML(3) + 1;
ETAsrc   = sqrt(URzz(nx_src, ny_src)/ERxx(nx_src,ny_src));
nsrc     = sqrt(URzz(nx_src, ny_src)*ERxx(nx_src,ny_src));
tau      = 0.5/f0;
t0       = nsrc*sqrt((Sx/2)^2 + (Sy/2)^2)/c0;

% CALCULATE NUMBER OF TIME STEPS
tprop = nmax*Sy/c0;
t     = 2*t0 + 2*tprop;
STEPS = ceil(t/dt);

% DIRECTION COSINES
gx = cos(theta);
gy = sin(theta);

% COMPENSATE FOR DISPERSION
k0     = 2*pi*f0/c0;
f      = c0*dt/sin(pi*f0*dt)*sqrt( (sin(k0*gx*dx/2)/dx)^2 ...
                                 + (sin(k0*gy*dy/2)/dy)^2);
                              
ERxx = f*ERxx;
ERyy = f*ERyy;
URzz = f*URzz;

% Q
[Y, X]  = meshgrid(ya-mean(ya), xa-mean(xa));
r       = 0.45*(Sx - NPML(1)*dx - NPML(2)*dx);
Q       = (X.^2 + Y.^2) >= r^2;

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
is2t(1:Nx-1,:) = Q(1:Nx-1,:) > Q(2:Nx,:);
it2s(1:Nx-1,:) = Q(1:Nx-1,:) < Q(2:Nx,:);
js2t(:,1:Ny-1) = Q(:,1:Ny-1) > Q(:,2:Ny);
jt2s(:,1:Ny-1) = Q(:,1:Ny-1) < Q(:,2:Ny);
isft(2:Nx,:)   = Q(2:Nx,:) > Q(1:Nx-1, :);
itfs(2:Nx,:)   = Q(2:Nx,:) < Q(1:Nx-1, :);
jsft(:,2:Ny)   = Q(:,2:Ny) > Q(:,1:Ny-1);
jtfs(:,2:Ny)   = Q(:,2:Ny) < Q(:,1:Ny-1);

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
sigx2 = zeros(Nx2, Ny2);
sigy2 = zeros(Nx2, Ny2);

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
    sigy2(:, ny1) = sigmax*(ny/2/NPML(3))^pml_Npml;
end

% ADD YHI CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(4)*dy2);
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy2(:, ny1) = sigmax*(ny/2/NPML(4))^pml_Npml;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M
M = c0*dt;

% CALCULATE UPDATE COEFFICIENTS Dx
sigx = sigx2(2:2:Nx2, 1:2:Ny2);
sigy = sigy2(2:2:Nx2, 1:2:Ny2);
bDxy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cDxy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bDxy - 1);

% CALCULATE UPDATE COEFFICIENTS Dy
sigx = sigx2(1:2:Nx2, 2:2:Ny2);
sigy = sigy2(1:2:Nx2, 2:2:Ny2);
bDyx = exp(-(sigx/pml_ky + pml_ax)*dt/e0);
cDyx = sigx./(sigx*pml_kx + pml_ax*pml_kx^2).*(bDyx - 1);

% CALCULATE UPDATE COEFFICIENTS Bz
sigx = sigx2(2:2:Nx2, 2:2:Ny2);
sigy = sigy2(2:2:Nx2, 2:2:Ny2);
bBzx = exp(-(sigx/pml_kx + pml_ax)*dt/e0);
cBzx = sigx./(sigx*pml_kx + pml_ax*pml_kx^2).*(bBzx - 1);
bBzy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cBzy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bBzy - 1);

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

% INITIALIZE CONVOLUTIONS

psiDx_ylo = zeros(Nx, NPML(3));
psiDx_yhi = zeros(Nx, NPML(4));

psiDy_xlo = zeros(NPML(1), Ny);
psiDy_xhi = zeros(NPML(2), Ny);

psiBz_xlo = zeros(NPML(1), Ny);
psiBz_xhi = zeros(NPML(2), Ny);
psiBz_ylo = zeros(Nx, NPML(3));
psiBz_yhi = zeros(Nx, NPML(4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN FDTD LOOP FOR H MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATE OVER TIME
%

for T = 1 : STEPS

    % CALCULATE SOURCE FIELDS
    Hzsrc =            exp(-(((T + 0.5)*dt - t0 - (gx* (X + dx/2) + gy*(Y + dy/2))*nsrc/c0)/tau).^2);
    Exsrc = -gy*ETAsrc*exp(-(((T + 0  )*dt - t0 - (gx* (X + dx/2) + gy*(Y       ))*nsrc/c0)/tau).^2);
    Eysrc = +gx*ETAsrc*exp(-(((T + 0  )*dt - t0 - (gx* (X       ) + gy*(Y + dy/2))*nsrc/c0)/tau).^2);
    
    % CALCULATE DERIVATIVES OF E
        % dEyx
        for ny = 1 : Ny
            for nx = 1 : Nx-1
                dEyx(nx, ny) = ( Ey(nx+1, ny) - Ey(nx,ny) )/dx;
            end
            dEyx(Nx, ny) = ( Ey(1, ny) - Ey(Nx,ny) )/dx;
        end
        
        % dExy
        for nx = 1 : Nx
            for ny = 1 : Ny - 1
                dExy(nx, ny) = ( Ex(nx, ny+1) - Ex(nx,ny) ) /dy; 
            end
            dExy(nx, Ny) = ( Ex(nx, 1) - Ex(nx,Ny) ) /dy; 
        end
     
     % INCORPORATE TF/SF CORRECTIONS
     dEyx(it2s) = dEyx(it2s) + Eysrc(it2s+1)/dx;
     dEyx(is2t) = dEyx(is2t) - Eysrc(is2t+1)/dx;
     dExy(jt2s) = dExy(jt2s) + Exsrc(jt2s+Nx)/dy;
     dExy(js2t) = dExy(js2t) - Exsrc(js2t+Nx)/dy;
     
     % UPDATE CONVOLUTIONS
     psiBz_xlo = bBzx(1:NPML(1),:).*psiBz_xlo ...
               + cBzx(1:NPML(1),:).*dEyx(1:NPML(1),:);
     psiBz_xhi = bBzx(Nx - NPML(2)+1:Nx, :).*psiBz_xhi ...
               + cBzx(Nx - NPML(2)+1:Nx, :).*dEyx(Nx - NPML(2)+1:Nx, :);
     psiBz_ylo = bBzy(:, 1:NPML(3)).*psiBz_ylo ...
               + cBzy(:, 1:NPML(3)).*dExy(:, 1:NPML(3)); 
     psiBz_yhi = bBzy(:, Ny - NPML(4)+1:Ny).*psiBz_yhi ...
               + cBzy(:, Ny - NPML(4)+1:Ny).*dExy(:, Ny - NPML(4)+1:Ny); 
     
           
     % UPDATE B FROM E
     Bz                         = Bz - M*(dEyx/pml_kx - dExy/pml_ky);
     Bz(1:NPML(1),:)            = Bz(1:NPML(1),:) - M*(psiBz_xlo);
     Bz(Nx - NPML(2)+1:Nx, :)   = Bz(Nx - NPML(2)+1:Nx, :) - M*(psiBz_xhi);
     Bz(:, 1:NPML(3))           = Bz(:, 1:NPML(3)) - M*(-psiBz_ylo);
     Bz(:,Ny-NPML(4)+1:Ny)      = Bz(:,Ny-NPML(4)+1:Ny) - M*(-psiBz_yhi);
     
     % UPDATE H FROM B
     Hz = Bz./URzz;
     
     % CALCULATE DERIVATIVES OF H
        % dHzy
        for nx = 1 : Nx
            dHzy(nx, 1) = ( Hz(nx,1) - Hz(nx, Ny) )/dy;
            for ny = 2 : Ny
                dHzy(nx, ny) = ( Hz(nx,ny) - Hz(nx, ny-1) )/dy;
            end
        end
        
        % dHzx
        for ny = 1 : Ny
            dHzx(1, ny) = ( Hz(1,ny) - Hz(Nx, ny) ) / dx;
            for nx = 2 : Nx
                dHzx(nx, ny) = ( Hz(nx,ny) - Hz(nx-1, ny) )/ dx;
            end
        end
      
     % INCORPORATE TF/SF CORRECTIONS
     dHzx(isft) = dHzx(isft) + Hzsrc(isft-1)/dx;
     dHzx(itfs) = dHzx(itfs) - Hzsrc(itfs-1)/dx;
     dHzy(jsft) = dHzy(jsft) + Hzsrc(jsft-Nx)/dy;
     dHzy(jtfs) = dHzy(jtfs) - Hzsrc(jtfs-Nx)/dy;
     
     % UPDATE CONVOLUTIONS
     psiDx_ylo = bDxy(:,1:NPML(3)).*psiDx_ylo ...
               + cDxy(:,1:NPML(3)).*dHzy(:,1:NPML(3));
     psiDx_yhi = bDxy(:, Ny - NPML(4)+1:Ny).*psiDx_yhi...
               + cDxy(:, Ny - NPML(4)+1:Ny).*dHzy(:, Ny - NPML(4)+1:Ny);
     psiDy_xlo = bDyx(1:NPML(1),:).*psiDy_xlo ... 
               + cDyx(1:NPML(1),:).*dHzx(1:NPML(1),:);      
     psiDy_xhi = bDyx(Nx-NPML(2)+1:Nx,:).*psiDy_xhi ... 
               + cDyx(Nx-NPML(2)+1:Nx,:).*dHzx(Nx-NPML(2)+1:Nx,:);      
     
        
     % UPDATE D FROM H
     Dx                         = Dx + M*(dHzy/pml_ky);
     Dx(:, 1:NPML(3))           = Dx(:, 1:NPML(3)) + M*(psiDx_ylo);
     Dx(:, Ny - NPML(4)+1:Ny)   = Dx(:, Ny - NPML(4)+1:Ny) + M*(psiDx_yhi);
     
     Dy                         = Dy + M*(-dHzx/pml_kx);
     Dy(1:NPML(1),:)            = Dy(1:NPML(1),:) + M*(-psiDy_xlo);
     Dy(Nx-NPML(2)+1:Nx,:)      = Dy(Nx-NPML(2)+1:Nx,:) + M*(-psiDy_xhi);
     
     
     % UPDATE E FROM D
     Ex = Dx./ERxx;
     Ey = Dy./ERyy;
     
     % SHOW FIELDS
     if mod(T, 10) == 0
         
         subplot(1,3,1);
         imagesc(xa, ya, Hz.')
         axis equal tight
         colorbar;
         title('Hz')
         caxis([-1 +1]);
         
         subplot(1,3,2);
         imagesc(xa, ya, Ex.')
         axis equal tight
         colorbar;
         title('Ex')
         caxis([-1 +1]);
         
         subplot(1,3,3);
         imagesc(xa, ya, Ey.')
         axis equal tight
         colorbar;
         title('Ey')
         caxis([-1 +1]);
         
         drawnow;
     end
end
















