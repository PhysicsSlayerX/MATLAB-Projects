% session2.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees = pi/180;

% OPEN FIGURE WINDOW
figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS OF THE PLANAR GRATING

a = 1;      %period of grating

% GRID PARAMETERS
Sx = 10*a;
Sy = 10*a;

NRESLO = 8;
NRESHI = 40;

theta = 45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOW RESOLUTION GRID

dx = a/NRESLO;
dy = a/NRESLO;

Nx = ceil(Sx/dx); 
dx = Sx/Nx;

Ny = ceil(Sy/dy);
dy = Sy/Ny;

xa = [0: Nx - 1]*dx; xa = xa - mean(xa);
ya = [0: Ny - 1]*dy; ya = ya - mean(ya);

[Y, X] = meshgrid(ya, xa);

% HIGH RESOLUTION GRID
dx2 = a/NRESHI;
dy2 = a/NRESHI;

Nx2 = ceil(Sx/dx2);
Ny2 = ceil(Sy/dy2);

xa2 = linspace(xa(1), xa(Nx), Nx2);
dx2 = xa2(2) - xa2(1);

ya2 = linspace(ya(1), ya(Ny), Ny2);
dy2 = ya2(2) - ya2(1);

[Y2, X2] = meshgrid(ya2, xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERIOD
PER = X - Y;
PER = PER - min(PER(:));
PER = PER / max(PER(:));
PER = 0.5*a + a*PER;

% ORIENTATION
r = 0.25*Sx;
RSQ = X.^2 + Y.^2;
THETA = (theta*degrees)*(RSQ < r^2);
THETA = svlblur(THETA, [2, 2]);

% FILL FACTOR
FF = X2 + Y2;
FF = FF - min(FF(:));
FF = FF / max(FF(:));
FF = 0.2 + 0.6*FF;


% SHOW PERIOD
subplot(2,3,1);
h = imagesc(xa, ya, PER');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
axis equal tight;
title('PER');

% SHOW PERIOD
subplot(2,3,2);
h = imagesc(xa, ya, THETA');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
axis equal tight;
title('THETA');

% SHOW FILL FRACTION
subplot(2,3,3);
h = imagesc(xa, ya, FF');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
axis equal tight;
title('Fill Fraction');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE SPATIALLY-VARIANT GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTRUCT K-FUNCTION
Kx = (2*pi./PER) .* cos(THETA);
Ky = (2*pi./PER) .* sin(THETA);

% Or you can use the svlsolve() function
PHI = svlsolve(Kx, Ky, dx, dy);

% INTERPOLATE TO HIGH RESOLUTION GRID
PHI = interp2(ya, xa', PHI, ya2, xa2');

% CALCULATE ANALOG GRATING
GA = cos(PHI);

% CALCULATE BINARY GRATING
GTH = cos(pi*FF);
GB = (GA > GTH);

% SHOW GRATING PHASE
subplot(2,3,4);
pcolor(xa2, ya2, PHI');
shading interp;
colorbar;
axis equal tight;
title('PHI');

% SHOW ANALOG GRATING
subplot(2,3,5);
pcolor(xa2, ya2, GA');
shading interp;
colorbar;
axis equal tight;
title('ANALOG GRATING');

% SHOW BINARY GRATING
subplot(2,3,6);
h = imagesc(xa2, ya2, GB');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
axis equal tight;
title('BINARY GRATING');








