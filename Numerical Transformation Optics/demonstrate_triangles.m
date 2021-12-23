% demonstrate_triangles.m

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


% TRIANGLE PARAMETERS
w = 9;

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
h = w*sqrt(3)/2;
ny = 1*round(h/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;


% INCORPORATE TRIANGLE
for ny = ny1 : ny2
    
    % Calculate the 0 to 1 number
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    
    % Calculate Width of Triangle at Current Position
    nx = round(f*w/dx);
    
    % Calculate Stop/Start Indices
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    
    % Incorporate 1's into A
    A(nx1:nx2,ny) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY ARRAY A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW A
imagesc(xa, ya, A.');

% SET VIEW
axis equal tight;
colorbar;





















