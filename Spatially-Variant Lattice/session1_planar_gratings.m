% session1.m

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

% PARAMETERS FOR THE PLANAR GRATING
a       = 0.4;               % period of grating
theta   = 120 * degrees;    % slant of the grating
ff      = 0.75;             % fill fraction

% GRID PARAMETERS
Sx      = 10;
Sy      = 10;
NRES    = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OUR GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = a/NRES;
dy = a/NRES;

% NUMBER OF POINTS

Nx = ceil(Sx/dx);
dx = Sx/Nx;

Ny = ceil(Sy/dy);
dy = Sy/Ny;

% GRID AXES
xa = [0: Nx - 1]*dx; xa = xa - mean(xa);
ya = [0: Ny - 1]*dy; ya = ya - mean(ya);

% MESHGRID
[Y, X] = meshgrid(ya, xa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE UNIFORM PLANAR GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRATING VECTOR FUNCTION
Kx = (2*pi/a) * cos(theta);
Ky = (2*pi/a) * sin(theta);

% CALCULATE ANALOG GRATING
GA = cos(Kx*X + Ky*Y);

% CALCULATE BINARY GRATING
gth = cos(pi*ff);
GB  = double(GA > gth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHOW GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW ANALOG GRATING
subplot(1, 2, 1);
pcolor(xa, ya, GA');
shading interp;
axis equal tight;
colorbar;
title('ANALOG GRATING');

% SHOW BINARY GRATING
subplot(1, 2, 2);
pcolor(xa, ya, GB');
shading interp;
axis equal tight;
colorbar;
title('BINARY GRATING');






