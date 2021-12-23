% demonstrate_rectangles.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID PARAMETERS
Sx = 10;
Sy = 10;

Nx = 100;
Ny = round(Nx*Sy/Sx);


% RECTANGLE PARAMETERS
wx = 7;
wy = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE A GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
dx = Sx/Nx;
dy = Sy/Ny;

% CALCULATE AXIS VECTORS
xa = [0.5:Nx-0.5]*dx;
ya = [0.5:Ny-0.5]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD RECTANGLE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE ARRAY TO ALL ZEROS
A = zeros(Nx, Ny);

% STOP/START INDICES
nx = round(wx/dx);
nx1 = 1 + floor((Nx - nx)/2);
nx2 = nx1 + nx - 1;

ny = round(wy/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;


% INCORPORATE RECTANGLE
A(nx1:nx2, ny1:ny2) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY ARRAY A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW A
imagesc(xa, ya, A.');

% SET VIEW
axis equal tight;
colorbar;





















