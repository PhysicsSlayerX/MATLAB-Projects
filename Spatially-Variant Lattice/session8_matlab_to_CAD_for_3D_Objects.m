% session8.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% START TIMER
time1 = clock; 

% UNITS
degrees = pi/180;

% OPEN FIGURE WINDOW
fig = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0.05 1 0.95]);
set(fig, 'Name', 'Generating Spatially-Variant Lattices');
set(fig, 'NumberTitle', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MOVIE CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Movie Creation
MAKE_MOVIE = 0;
movie_title = 'Generating_3D_SV_Lattice.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    vidObj.FrameRate = 30;
    VidObj.Quality = 100;
    open(vidObj);
end

% PLOT CONTROL
twoD_Slice = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD DATA FROM FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% READ BINARY LATTICE FROM FILE
load svlattice3D SVLB dx2 dy2 dz2;

% CALCULATE GRID
[Nx, Ny, Nz] = size(SVLB);
dx = dx2;
dy = dy2;
dz = dz2;

Sx = Nx*dx;
Sy = Ny*dy;
Sz = Nz*dz;

xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;
za = [0:Nz-1]*dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MESH LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORCE TOP AND BOTTOM OF LATTICE TO BE THE SAME
S = ( SVLB (:,:,1) + SVLB(:,:,Nz))/2;
SVLB(:,:,1)  = S;
SVLB(:,:,Nz) = S;

% SMOOTH LATTICE
SVL = smooth3(SVLB);

if twoD_Slice
    % 2D SLICE
    subplot(1,2,1);
    nzc = round(Nz/2);
    imagesc(xa, ya, SVLB(:,:, nzc)');
    axis equal tight;
    colorbar;

    subplot(1,2,2);
    nzc = round(Nz/2);
    imagesc(xa, ya, SVL(:,:, nzc)');
    axis equal tight;
    colorbar;
end

% MESH THE SURFACE
[F, V]    = isosurface(ya, xa, za, SVL, 0.5);
[FC,  VC] = isocaps(ya, xa, za, SVL, 0.5);

% COMBINE FACES AND VERTICES
F = [F ; FC+length(V(:,1))];
V = [V ; VC];
clear FC VC;

% DRAW MESH
h = patch('faces', F, 'vertices', V, 'FaceColor', 'g', 'LineStyle', 'none');
view(10, 60);
axis equal tight;
camlight; lighting phong;
title('SPATIALLY-VARIANT LATTICE');
drawnow;

% SAVE LATTICE TO STL FILE
disp('Creating STL File...'); drawnow;
svlcad('svlattice3D.stl', F, V);


if MAKE_MOVIE == 1
    close(vidObj);
end


% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);






