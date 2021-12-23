% FDTD2D_HMode_Step6_DiffractionGrating.m
% Simulate scattering from a sawtooth diffraction grating

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
f0      = 10 * gigahertz;
lam0    = c0/f0;
k0      = 2*pi/lam0;
bw      = 5*lam0;
theta   = 0*degrees;

% DIFFRACTION GRATING PARAMETERS
L = 1.5 * centimeters;
d = 1.0 * centimeters;
er = 9.0;
NP = 50;


% PML PARAMETERS
pml_kx      = 1;
pml_ky      = 1;
pml_ax      = 1e-10;
pml_ay      = 1e-10;
pml_Npml    = 3;
pml_R0      = 1e-8;

% GRID PARAMETERS
NRES     = 20;
Sy       = 20 * lam0;
NPML     = [ 20 20 20 20 ];
nmax     = sqrt(er);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);
dx = L/nx;
ny = ceil(d/dy);
dy = d/ny;

% COMPUTER GRID SIZE
Sx = NP*L;
Nx = NPML(1) + ceil(Sx/dx) + NPML(2);
Sx = Nx*dx;
Ny = NPML(3) + ceil(Sy/dy) + NPML(4);
Sy = Ny*dy;

% GRID AXES
xa = [0 : Nx - 1]*dx;
ya = [0 : Ny - 1]*dy;
[Y, X] = meshgrid(ya, xa);

% 2X GRID
Nx2 = 2 * Nx; dx2 = dx/2;
Ny2 = 2 * Ny; dy2 = dy/2;
xa2 = [0:Nx2-1]*dx2;
xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2;
ya2 = ya2 - mean(ya2);
[Y2, X2] = meshgrid(ya2,xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO VACUUM
ER2 = ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);

% BUILD ONE UNIT CELL
ny = round(d/dy2);
ny1 = 1 + floor((Ny2 - ny)/2);
ny2 = ny1 + ny - 1;
for ny = ny1 : ny2
    f   = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx  = round(f*L/dx2);
    nx2 = round(L/dx2);
    nx1 = nx2 - nx + 1;
    ER2(nx1:nx2, ny) = er;
end

% STACK LEFT TO RIGHT
nx = round(L/dx2);
nxa1 = 1;
nxa2 = nxa1 + nx - 1;
nxb1 = nxa2 + 1;
nxb2 = nxb1 + nx - 1;
while nxb2 <= Nx2
    ER2(nxb1:nxb2, ny1:ny2) = ER2(nxa1:nxa2, ny1:ny2);
    nxb1 = nxb2 + 1;
    nxb2 = nxb1 + nx - 1;
end

nx = Nx2 - nxb1 + 1;
ER2(nxb1:Nx2, ny1:ny2) = ER2(nxa1:nxa1+nx-1,ny1:ny2);

% ADD SUBSTRATE
ER2(:, ny2+1:Ny2) = er;

% EXTRACT MATERIALS ARRAYS FROM 2X GRID
ERxx = ER2(2:2:Nx2, 1:2:Ny2);
ERyy = ER2(1:2:Nx2, 2:2:Ny2);
URzz = UR2(2:2:Nx2, 2:2:Ny2);

% % SHOW DEVICE
% imagesc(xa,ya2,ERxx.');
% axis equal tight
% colorbar;
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE STABLE TIME STEP
dmin     = min([dx dy]);
dt       = dmin/(2*c0);

% CALCULATE MATERIALS WHERE SOURCE IS INJECTED
nx_src   = 1;
ny_src   = 1;
ETAsrc   = sqrt(URzz(nx_src, ny_src)/ERxx(nx_src,ny_src));
nsrc     = sqrt(URzz(nx_src, ny_src)*ERxx(nx_src,ny_src));

% CALCULATE NUMBER OF TIME STEPS
tprop = nmax*Sy/c0;
t     = 2*tprop;
STEPS = ceil(t/dt);

% CALCULATE FREQUENCY-DOMAIN GAUSSIAN BEAM
[TH, R] = cart2pol(X2, Y2);
[XR, YR] = pol2cart(TH+theta, R);
G_beam = exp(-(2*XR/bw).^2);
G_wave = exp(-1i*k0*nsrc*YR);
G = G_beam.*G_wave .* (Y2 <= 0);

% CALCULATE RAMP FUNCTION
t = [0:STEPS-1]*dt;
dur = 10/f0;
RMP = (1 + sin(0.5*pi/dur*(2*t-dur)))/2;
RMP(t >= dur) = 1;
% 
% plot(t, RMP);
% return

% imagesc(xa2, ya2, real(G).');
% axis equal tight;
% colorbar;
% return

% COMPENSATE FOR DISPERSION
gx = sin(theta);
gy = cos(theta);
k0     = 2*pi*f0/c0;
f      = c0*dt/sin(pi*f0*dt)*sqrt( (sin(k0*gx*dx/2)/dx)^2 ...
                                 + (sin(k0*gy*dy/2)/dy)^2);
                              
ERxx = f*ERxx;
ERyy = f*ERyy;
URzz = f*URzz;

% Q

nx1 = NPML(1) + 2;
nx2 = Nx - NPML(2) - 1;
ny1 = NPML(3) + 2;
ny2 = Ny - NPML(4) - 1;
Q   = ones(Nx,Ny);
Q(nx1:nx2,ny1:ny2) = 0;

% [Y, X]  = meshgrid(ya-mean(ya), xa-mean(xa));
% r       = 0.45*(Sx - NPML(1)*dx - NPML(2)*dx);
% Q       = (X.^2 + Y.^2) >= r^2;



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
    Hzsrc =            RMP(T)*real(G(2:2:Nx2, 2:2:Ny2)*exp(1i*2*pi*f0*(T + 0.5)*dt));
    Exsrc = -gy*ETAsrc*RMP(T)*real(G(2:2:Nx2, 1:2:Ny2)*exp(1i*2*pi*f0*(T      )*dt));
    Eysrc = +gx*ETAsrc*RMP(T)*real(G(1:2:Nx2, 2:2:Ny2)*exp(1i*2*pi*f0*(T      )*dt));
    
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
     if mod(T, 100) == 0
         
         %Refresh Graphics
         clf;
         hold on;
         
         
         % Show Field
         imagesc(xa, ya, Hz.');
         caxis([-1 +1]);
         
         % Add Device
         contour(X,Y, ERxx, 8.5, 'Color','w');
         
         % Set View
         hold off;
         axis equal tight off;
         set(gca, 'YDir', 'reverse');
         colorbar;
         
         
         drawnow;
     end
end
















