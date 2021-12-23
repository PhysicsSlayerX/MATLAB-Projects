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
NFREQ   = 100;
freq1   = 1* gigahertz;
freq2   = 5 * gigahertz;
FREQ    = linspace(freq1, freq2, NFREQ);

% LOSSY SLAB PARAMETERS

f0      = mean(FREQ);
lam0    = c0/f0;
% DIFFRACTION GRATING PARAMETERS
L       = 1.5 * centimeters;
d       = 1.0 * centimeters;
er      = 15;
sig     = 100;

% PML PARAMETERS
pml_ky   = 1;
pml_ay   = 1e-10;
pml_Npml = 3;
pml_R0   = 1e-8;

% GRID PARAMETERS
lam_max = c0/min(FREQ);
NRES    = 28;
SPACER  = 0.05*lam_max * [1 1];
NPML    = [20 20];
nmax    = sqrt(er);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% SET FIGURE WINDOW
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Color','k', 'Position', [0 42 1922 954]);
set(fig, 'Name', 'FDTD-2D Frequency Response Graph');
set(fig, 'NumberTitle', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
lam_min     = c0/max(FREQ);
dx          = lam_min/nmax/NRES;
dy          = lam_min/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);
dx = L/nx;
ny = ceil(d/dy);
dy = d/ny;

% COMPUTE GRID SIZE
Sx = L;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + d + SPACER(2);
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

ERzz    = ones(Nx,Ny);
URxx    = ones(Nx,Ny);
URyy    = ones(Nx,Ny);
SIGzz   = zeros(Nx,Ny);

addpath 'D:\Research\Matlab practice\Numerical Transformation Optics\SNN_images'
imagefiles = dir('D:\Research\Matlab practice\Numerical Transformation Optics\SNN_images\*.png');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   images{ii} = currentimage;
end

%nfiles = 5;
for imageIndex = 1:nfiles
    %imageIndex = 400;
    thisImage = imrotate(images{imageIndex}, 90, 'bilinear', 'crop');
    thisImage = imresize(thisImage,[Nx Nx]);

    lastimage = ceil(thisImage(:,:,1));
    lastimage = im2double(lastimage);

    % BUILD IMAGE ARRAY
    ny1 = NPML(1) + round(0.75*SPACER(1)/dy) + 1;
    ny2 = ny1 + round(L/dx) - 1;

    % ADD PERMITTIVITY

    lastimage(lastimage > 0) = ceil(lastimage(lastimage > 0));
    ERzz(:,ny1:ny2) = er*lastimage;
    SIGzz(:,ny1:ny2) = sig*lastimage;

    ERzz = flip(ERzz, 1);
    SIGzz = flip(SIGzz, 1);
    ERzz(ERzz <= 1) = 1;
    SIGzz(SIGzz <= 0) = 0;

    %ERzz(ERzz > 1) = er;

    % SHOW DEVICE ERzz
    subplot(1, 2, 1);
    imagesc(xa, ya, ERzz.');
    axis equal tight off;
    title('ERzz');
    colorbar;
    plot_darkmode

    % SHOW DEVICE SIGzz
    subplot(1, 2, 2);
    imagesc(xa, ya, SIGzz.');
    axis equal tight off;
    title('SIGzz');
    colorbar;
    plot_darkmode;
    
    %drawnow;
    % EXTRACT ER and SIG MATRICES 
    imageER_zz      = mat2gray(ERzz(:,ny1:ny2));
    imageSIG_zz     = mat2gray(SIGzz(:,ny1:ny2));
end























