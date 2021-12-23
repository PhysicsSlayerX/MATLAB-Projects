% session5.m

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
movie_title = 'Generating_SV_Lattice.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    vidObj.FrameRate = 30;
    VidObj.Quality = 100;
    open(vidObj);
end

% PLOT CONTROL
plotControl = 1;
plot_now = 0;

% SAVE STL FILE
save_stl = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD DATA FROM FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD BINARY LATTICE FROM FILE
load svlattice2D  SVLB dx2 dy2;
SVLB = 1 - SVLB;

% RECALCULATE GRID
[Nx, Ny] = size(SVLB);
dx = dx2;
dy = dy2;

Sx = Nx*dx;
Sy = Ny*dy;

xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;

% SHOW BINARY LATTICE
subplot(1,3,1);
hh = imagesc(xa, ya, SVLB');
h2 = get(hh, 'Parent');
set(h2, 'YDir' , 'normal');
axis equal tight;
title('BINARY LATTICE');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REDUCE RESOLUTION OF LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOWER RESOLUTION GRID
Nx2 = 128;
Ny2 = round(Nx2*Sy/Sx);

xa2 = linspace(xa(1), xa(Nx), Nx2);
ya2 = linspace(ya(1), ya(Ny), Ny2);

dx2 = xa2(2) - xa2(1);
dy2 = ya2(2) - ya2(1);

% INTERPOLATE TO LOWER RESOLUTION GRID
SVL = svlblur(SVLB, [Nx, Ny]./[Nx2 Ny2]);
SVL = interp2(ya, xa', SVL, ya2, xa2');

% SHOW LOWER RESOLUTION LATTICE
subplot(1,3,2);
hh = imagesc(xa, ya, SVL');
h2 = get(hh, 'Parent');
set(h2, 'YDir' , 'normal');
axis equal tight;
title('REDUCED LATTICE');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MESH THE LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STACK SVL
SVL = rot90(fliplr(SVL));
SVL(:, :, 2) = SVL;

% GENERATE MESH
[F, V] = isocaps(xa2, ya2,[0 1], SVL, 0.5, 'zmax');

% SHOW MESH
subplot(1,3,3);
c = [0.6 0.6 1.0];
h = patch('faces', F, 'vertices', V, 'FaceColor', c);
axis equal tight;
view(0, 90);
title('MESH');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE TO STL FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE STL FILE
if save_stl
    svlcad('svlattice2D.stl', F, V);
end

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






