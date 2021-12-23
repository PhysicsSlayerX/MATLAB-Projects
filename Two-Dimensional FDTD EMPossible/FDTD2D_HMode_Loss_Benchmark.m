% FDTD2D_HMode_Loss_Benchmark.m
% Implement Loss

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
NFREQ = 100;
freq1 = 5*gigahertz;
freq2 = 15*gigahertz;
FREQ  = linspace(freq1,freq2, NFREQ);

% LOSSY SLAB PARAMETERS
f0      = mean(FREQ);
lam0    = c0/f0;
d       = 0.25*lam0;
er      = 9.0;
sig     = 100;

% PML PARAMETERS
pml_ky      = 1;
pml_ay      = 1e-10;
pml_Npml    = 3;
pml_R0      = 1e-8;

% GRID PARAMETERS
NRES     = 50;
Sx       = 1*lam0;
SPACER   = 0.5 * lam0 * [ 2 2 ];
NPML     = [ 20 20 ];
nmax     = sqrt(er);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
lam_min     = c0/max(FREQ);
dx          = lam_min/nmax/NRES;
dy          = lam_min/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
ny = ceil(d/dy);
dy = d/ny;


% COMPUTER GRID SIZE
Nx = ceil(Sx/dx);
Sx = Nx*dx;
Sy = SPACER(1) + d + SPACER(2);
Ny = NPML(1) + ceil(Sy/dy) + NPML(2);
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
ER2 = ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);
SIG2 = zeros(Nx2, Ny2);

% BUILD SLAB
ny1 = 2*NPML(1) + round(SPACER(1)/dy2) + 1;
ny2 = ny1 + round(d/dy2) - 1;
ER2(:,ny1:ny2) = er;
SIG2(:, ny1:ny2) = sig;

% EXTRACT MATERIAL ARRAYS FROM 2X GRID
ERxx  = ER2(2:2:Nx2, 1:2:Ny2);
ERyy  = ER2(1:2:Nx2, 2:2:Ny2);
URzz  = UR2(2:2:Nx2, 2:2:Ny2);
SIGxx = SIG2(2:2:Nx2,1:2:Ny2);
SIGyy = SIG2(1:2:Nx2,2:2:Ny2);

% SHOW DEVICE
subplot(1, 5, 1);
imagesc(xa,ya, ERxx.');
axis equal tight off;
title('ERxx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE STABLE TIME STEP
dmin     = min([dx dy]);
dt       = dmin/(2*c0);

% POSITION OF SOURCE
ny_ref = NPML(1) + 1;
ny_src = ny_ref + 1;
ny_trn = Ny - NPML(2);

% MATERIAL PROPERTIES
ersrc = ERxx(1, ny_ref);
ertrn = ERxx(1, ny_trn);
nref  = sqrt(URzz(1, ny_ref)*ersrc);
ntrn  = sqrt(URzz(1, ny_ref)*ertrn);

% CALCULATE GAUSSIAN PULSE PARAMETERS
tau      = 0.5/max(FREQ);
t0       = 3*tau;

% CALCULATE NUMBER OF TIME STEPS
tprop = nmax*Sy/c0;
t     = 2*t0 + 15*tprop;
STEPS = ceil(t/dt);

% CALCULATE GAUSSIAN SOURCE
ETAsrc   = sqrt(URzz(1, ny_src)/ERxx(1,ny_src));
nsrc     = sqrt(URzz(1, ny_src)*ERxx(1,ny_src));
delt     = 0.5*dt + 0.5*nsrc*dy/c0;
t        = [0:STEPS-1]*dt;
Hzsrc    =         exp(-((t - t0       )/tau).^2);
Exsrc    = -ETAsrc*exp(-((t - t0 - delt)/tau).^2);

% COMPENSATE FOR DISPERSION
f0     = mean(FREQ);
k0     = 2*pi*f0/c0;
f      = c0*dt/sin(pi*f0*dt)*sin(k0*dy/2)/dy;
                              
ERxx = f*ERxx;
ERyy = f*ERyy;
URzz = f*URzz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PML CONDUCTIVITY ON 2X GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE PML CONDUCTIVITIES
sigy2 = zeros(Nx2, Ny2);

% ADD YLO CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(1)*dy2);
for ny = 1 : 2*NPML(1)
    ny1 = 2*NPML(1) - ny + 1;
    sigy2(:, ny1) = sigmax*(ny/2/NPML(1))^pml_Npml;
end

% ADD YHI CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(2)*dy2);
for ny = 1 : 2*NPML(2)
    ny1 = Ny2 - 2*NPML(2) + ny;
    sigy2(:, ny1) = sigmax*(ny/2/NPML(2))^pml_Npml;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M
M = c0*dt;

% CALCULATE UPDATE COEFFICIENTS Dx
sigy = sigy2(2:2:Nx2, 1:2:Ny2);
bDxy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cDxy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bDxy - 1);


% CALCULATE UPDATE COEFFICIENTS Bz
sigy = sigy2(2:2:Nx2, 2:2:Ny2);
bBzy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cBzy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bBzy - 1);

% CALCULATE UPDATE COEFFICIENTS FOR Ex and Ey
Ax = (c0*dt)*SIGxx;
Ay = (c0*dt)*SIGyy;

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

psiDx_ylo = zeros(Nx, NPML(1));
psiDx_yhi = zeros(Nx, NPML(2));
psiBz_ylo = zeros(Nx, NPML(1));
psiBz_yhi = zeros(Nx, NPML(2));

% INITIALIZE INTEGRATIONS
IEx = zeros(Nx, Ny);
IEy = zeros(Nx, Ny);

% INITIALIZE FOURIER TRANSFORMS
K     = exp(-1i*2*pi*FREQ*dt);
HREF  = zeros(Nx, NFREQ);
HTRN  = zeros(Nx, NFREQ);
SRC   = zeros(1, NFREQ);

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
     dExy(:,ny_src-1) = dExy(:,ny_src-1) - Exsrc(T)/dy;
     
     % UPDATE CONVOLUTIONS
     psiBz_ylo = bBzy(:, 1:NPML(1)).*psiBz_ylo ...
               + cBzy(:, 1:NPML(1)).*dExy(:, 1:NPML(1)); 
     psiBz_yhi = bBzy(:, Ny - NPML(2)+1:Ny).*psiBz_yhi ...
               + cBzy(:, Ny - NPML(2)+1:Ny).*dExy(:, Ny - NPML(2)+1:Ny); 
     
           
     % UPDATE B FROM E
     Bz                         = Bz - M*(dEyx - dExy/pml_ky);
     Bz(:, 1:NPML(1))           = Bz(:, 1:NPML(1)) - M*(-psiBz_ylo);
     Bz(:,Ny-NPML(2)+1:Ny)      = Bz(:,Ny-NPML(2)+1:Ny) - M*(-psiBz_yhi);
     
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
     dHzy(:, ny_src) = dHzy(:,ny_src) - Hzsrc(T)/dy;
     
     % UPDATE CONVOLUTIONS
     psiDx_ylo = bDxy(:,1:NPML(1)).*psiDx_ylo ...
               + cDxy(:,1:NPML(1)).*dHzy(:,1:NPML(1));
     psiDx_yhi = bDxy(:, Ny - NPML(2)+1:Ny).*psiDx_yhi...
               + cDxy(:, Ny - NPML(2)+1:Ny).*dHzy(:, Ny - NPML(2)+1:Ny);     
        
     % UPDATE D FROM H
     Dx                         = Dx + M*(dHzy/pml_ky);
     Dx(:, 1:NPML(1))           = Dx(:, 1:NPML(1)) + M*(psiDx_ylo);
     Dx(:, Ny - NPML(2)+1:Ny)   = Dx(:, Ny - NPML(2)+1:Ny) + M*(psiDx_yhi);

     Dy                         = Dy + M*(-dHzx);
     
     % UPDATE E FROM D
     IEx = IEx + Ex;
     IEy = IEy + Ey;
     
     Ex = (Dx - Ax.*IEx)./ERxx;
     Ey = (Dy - Ay.*IEy)./ERyy;
     
     % UPDATE FOURIER TRANSFORMS
     for nfreq = 1 : NFREQ
        HREF(:, nfreq) = HREF(:, nfreq) + (K(nfreq)^T)*Hz(:,ny_ref);
        HTRN(:, nfreq) = HTRN(:, nfreq) + (K(nfreq)^T)*Hz(:,ny_trn);
        SRC(nfreq)     = SRC(nfreq)     + (K(nfreq)^T)*Hzsrc(T);
     end
     
     % SHOW FIELDS
     if mod(T, 10) == 0
         
         % Calculate REF & TRN
         REF = zeros(1, NFREQ);
         TRN = zeros(1, NFREQ);
         for nfreq = 1 : NFREQ
            % get next frequency
            lam0    = c0/FREQ(nfreq);
            k0      = 2*pi*lam0;
            % wave vector components
            m       = [-floor(Nx/2):+floor((Nx-1)/2)].';
            kxi     = 0 - m*2*pi/Sx;
            kyinc   = k0*nref;
            kyref  = sqrt((k0*nref)^2 - kxi.^2);
            kytrn  = sqrt((k0*ntrn)^2 - kxi.^2);
            % reflection
            href = HREF(:, nfreq)/SRC(nfreq);
            aref = fftshift(fft(href))/Nx;
            RDE  = abs(aref).^2.*real(kyref/kyinc);
            REF(nfreq) = sum(RDE);
            % transmission
            htrn = HTRN(:,nfreq)/SRC(nfreq);
            atrn = fftshift(fft(htrn))/Nx;
            TDE  = abs(atrn).^2.*real(ersrc/ertrn*kytrn/kyinc);
            TRN(nfreq) = sum(TDE);
            
         end
         CON = REF + TRN;
         
         subplot(1,5,2);
         imagesc(xa, ya, Hz.')
         axis equal tight
         title('Hz')
         caxis([-1 +1]);
         
         % Show Frequency Response
         subplot(1,5,3:5);
         plot(FREQ/gigahertz, CON, 'k:');
         hold on;
         plot(FREQ/gigahertz, REF, '-r');
         plot(FREQ/gigahertz, TRN, '-b');
         hold off;
         axis tight;
         xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
         ylim([-0.05 1.1]);
         xlabel('Frequency (GHz)');
         title(['Iteration ' num2str(T) ' of ' num2str(STEPS)]);
         
         drawnow;
     end
end
















