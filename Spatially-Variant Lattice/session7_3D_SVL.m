% session7.m

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
plotControl = 1;
plot_now = 0;

% GENERATE FILL FRACTION SWEEP CONTROL
sweepControl = 0;
save_fill_fraction_data = 0;

% SHOW LATTICE FORMATION
show_lattice_formation = 0;

% SAVE BINARY LATTICE FILE
save_binary_lattice = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CEL DIMENSIONS;
a = 1;

% GRID PARAMETERS
Sx = 7*a;
Sy = 7*a;
Sz = 1*a;
NRESLO = 10;
NRESHI = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOW RESOLUTION GRID
dx = a/NRESLO;
dy = a/NRESLO;
dz = a/NRESLO;

Nx = ceil(Sx/dx); dx = Sx/Nx;
Ny = ceil(Sy/dy); dy = Sy/Ny;
Nz = ceil(Sz/dz); dz = Sz/Nz;

xa = [0: Nx - 1]*dx; xa = xa - mean(xa);
ya = [0: Ny - 1]*dy; ya = ya - mean(ya);
za = [0: Nz - 1]*dz; za = za - mean(za);

[Y, X, Z] = meshgrid(ya, xa, za);

% HIGH RESOLUTIOn GRID
dx2 = a/NRESHI;    Nx2 = ceil(Sx/dx2);
dy2 = a/NRESHI;    Ny2 = ceil(Sy/dy2);
dz2 = a/NRESHI;    Nz2 = ceil(Sz/dz2);

xa2 = linspace(xa(1), xa(Nx), Nx2);
ya2 = linspace(ya(1), ya(Ny), Ny2);
za2 = linspace(za(1), za(Nz), Nz2);

dx2 = xa2(2) - xa2(1);
dy2 = ya2(2) - ya2(1);
dz2 = za2(2) - za2(1);

[Y2, X2, Z2] = meshgrid(ya2,xa2,za2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE LIST OF GRATING VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIMPLE CUBIC UNIT CELL
K1 = (2*pi/a) * [1;0;0];
K2 = (2*pi/a) * [0;1;0];
K3 = (2*pi/a) * [0;0;1];

% GENERATE LIST
KLIST = [ K1 K2 K3 ];
CLIST = [  1  1  1 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FILL FRACTION SWEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL GRID
Nxu = ceil(a/dx2);          dxu = a/Nxu;
Nyu = ceil(a/dy2);          dyu = a/Nyu;
Nzu = ceil(a/dz2);          dzu = a/Nzu;

xau = [0:Nxu-1]*dxu;        xau = xau - mean(xau);
yau = [0:Nyu-1]*dyu;        yau = yau - mean(yau);
zau = [0:Nzu-1]*dzu;        zau = zau - mean(zau);

[YU, XU, ZU] = meshgrid(yau, xau, zau);

% BUILD THE ANALOG UNIT CELL
NK = length(CLIST);
UC = zeros(Nxu, Nyu, Nzu);
for nk = 1 : NK
    Kx = KLIST(1, nk);
    Ky = KLIST(2, nk);
    Kz = KLIST(3, nk);
    GA = exp(1i*(Kx*XU + Ky*YU + Kz*ZU));
    UC = UC + CLIST(nk)*GA;
end

% CLEAN UP NUMERICAL NOISE
UC = real(UC);
UC = UC - min(UC(:));
UC = UC / max(UC(:));

% GENERATE FILL FRACTION SWEEP

if sweepControl
    
    % Initialize Threshold Data
    gth_dat = linspace(0,1,20);
    
    % Tabulate Fill Fraction Data
    ff_dat = 0*gth_dat;
    for n = 1 : length(gth_dat)
        % Generate Binary Unit Cell
        UCB = UC > gth_dat(n);
        
        % calculate fill fraction
        ff_dat(n) = sum(UCB(:))/(Nxu*Nyu*Nzu);
        
        % Show Unit Cell
        clf;
        subplot(1,2,1);
        UCS = smooth3(UCB);
        s = isosurface(YU, XU, ZU, UCS, 0.5);
        h = patch(s, 'FaceColor', 'g');
        s = isocaps(YU, XU, ZU, UCS, 0.5);
        h = patch(s, 'FaceColor', 'g');
        view(10,10);
        axis equal tight;
        xlim([-a/2 a/2]);
        ylim([-a/2 a/2]);
        zlim([-a/2 a/2]);
        camlight; lighting phong;
        title('UNIT CELL');
        
        subplot(1,2,2);
        plot(gth_dat(1:n), ff_dat(1:n));
        xlim([gth_dat(1) max(gth_dat)]);
        ylim([0 1]);
        xlabel('Threshold Value');
        ylabel('Fill Fraction');
        
        drawnow;
    end
    
    if save_fill_fraction_data
        % Save to File
        save ffdat gth_dat ff_dat;
    end
% LOAD FILL FRACTION SWEEP FROM FILE
else
    
    % Load Data From File 
    load ffdat gth_dat ff_dat;
   
    % Generate Binary Unit Cell
    UCB = UC > 0.6;
   
    % Show Unit Cell
    clf;
    subplot(2,3,1);
    UCS = smooth3(UCB);
    s = isosurface(YU, XU, ZU, UCS, 0.5);
    h = patch(s, 'FaceColor', 'g');
    s = isocaps(YU, XU, ZU, UCS, 0.5);
    h = patch(s, 'FaceColor', 'g');
    view(10,10);
    axis equal tight;
    xlim([-a/2 a/2]);
    ylim([-a/2 a/2]);
    zlim([-a/2 a/2]);
    camlight; lighting phong;
    title('UNIT CELL');

    subplot(2,3,[2 3]);
    plot(gth_dat, ff_dat);
    xlim([gth_dat(1) max(gth_dat)]);
    ylim([0 1]);
    xlabel('Threshold Value');
    ylabel('Fill Fraction');

    drawnow; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERIOD
PER = a*ones(Nx, Ny, Nz);

% ORIENTATION (LOW RESOLUTION GRID)
Rz = atan2(Y + 0.51*Sy, X + 0.51*Sx); % Rotation About z-axis

% FILL FRACTION (HIGH RESOLUTION GRID)
FF = 0.35*ones(Nx2, Ny2, Nz2);

% SHOW INPUT DATA
nzc = round(Nz/2);

subplot(2,3,4);
h  = imagesc(xa, ya, PER(:, :, nzc)');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
title('PER');

subplot(2,3,5);
h  = imagesc(xa, ya, Rz(:, :, nzc)/degrees');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
title('RZ');

subplot(2,3,6);
h  = imagesc(xa, ya, FF(:, :, nzc)');
h2 = get(h, 'Parent');
set(h2, 'YDir', 'normal');
colorbar;
title('RZ');

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE SPATIALLY-VARIANT 3D LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE CENTER OF THE LOW RESOLUTION GRID
nxc = round(Nx/2);
nyc = round(Ny/2);
nzc = round(Nz/2);

%
% MAIN LOOP -- ITERATE OVER FOURIER EXPANSION
%

NK = length(CLIST);
SVL = zeros(Nx2, Ny2, Nz2);
for nk = 1 : NK
    
    % Construct K-Function
    KX = zeros(Nx, Ny, Nz);
    KY = zeros(Nx, Ny, Nz);
    KZ = zeros(Nx, Ny, Nz);
    for nz = 1 : Nz
        for ny = 1 : Ny
            for nx = 1 : Nx
                % Construct rotation matrix
                theta = Rz(nx, ny, nz);
                R = [cos(theta) -sin(theta) 0 ;...
                     sin(theta)  cos(theta) 0 ;...
                     0           0          1];
                 
                % Rotate K vector
                Ki = R*KLIST(:, nk);
                
                % record Ki
                KX(nx, ny, nz) = Ki(1);
                KY(nx, ny, nz) = Ki(2);
                KZ(nx, ny, nz) = Ki(3);
                 
            end
        end
    end
    
    % Calculate Grating Phase
    PHI = svlsolve(KX, KY, KZ, dx, dy, dz);
    PHI = PHI - PHI(nxc, nyc, nzc);
    
    % Interpolate to High Resolutin Grid
    PHI = interp3(ya, xa', za, PHI, ya2, xa2', za2);
    
    % Calcualte Analog Planar Grating
    GA = exp(1i*PHI);
    
    % Add Planar Grating to Overall Lattice
    SVL = SVL + CLIST(nk)*GA;
    
    % Show Graphical Status
    
    if show_lattice_formation
        clf;
        subplot(1,2,1);
        h = imagesc(xa2, ya2, real(GA(:,:, nzc)'));
        h2 = get(h, 'Parent');
        set(h2, 'YDir' , 'normal');
        axis equal tight;
        colorbar;
        title('CURRENT PLANAR GRATING');
        
        subplot(1,2,2);
        h = imagesc(xa2, ya2, real(SVL(:,:, nzc)'));
        h2 = get(h, 'Parent');
        set(h2, 'YDir' , 'normal');
        axis equal tight;
        colorbar;
        title('ANALOG LATTICE');
        
        drawnow;
    end
         
end

% CLEAN UP NUMERICAL NOISE
SVL = real(SVL);
SVL = SVL - min(SVL(:));
SVL = SVL / max(SVL(:));

% GENERATE BINARY LATTICE
GTH = interp1(ff_dat, gth_dat, FF(:));
GTH = reshape(GTH, Nx2, Ny2, Nz2);
SVLB = SVL > GTH;

% SHOW 3D LATTICE

clf;
SVLS = smooth3(SVLB);
s = isosurface(Y2, X2, Z2, SVLS, 0.5);
h = patch(s, 'FaceColor', 'g');
s = isocaps(Y2, X2, Z2, SVLS, 0.5);
h = patch(s, 'FaceColor', 'g');
view(10,60);
axis equal tight;
camlight; lighting phong;
title('SPATIALLY-VARIANT LATTICE');

% SAVE BINARY LATTICE TO FILE
if save_binary_lattice
    save svlattice3D SVLB dx2 dy2 dz2;
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






