% session3.m

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
plotLatticeFormation = 1;
plot_now = 0;

% SAVE LATTICE
save_lattice = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL PARAMETERS
a      = 1;

% GRID PARAMETERS
Nxh = 255;              % unit cell grid
Nyh = 255;

Sx = 10*a;
Sy = 10*a;
NRESLO = 10;            % number of points per period
NRESHI = 20;

% FOURIER EXPANSION PARAMETERS
cth = 0.005;
NP  = 51;
NQ  = 51;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOW RESOLUTION GRID
dx = a/NRESLO;
dy = a/NRESLO;

Nx = ceil(Sx/dx);   dx = Sx/Nx;
Ny = ceil(Sy/dy);   dy = Sy/Ny;

xa = [1:Nx]*dx;
ya = [1:Ny]*dy;

[Y, X] = meshgrid(ya, xa);

% HIGH RESOLUTION GRID
Kmax = (2*pi/a) * [floor(NP/2); floor(NQ/2)]; % max grating vector
amin = 2*pi/norm(Kmax);                       % minumum period
dx2  = amin/NRESHI;
dy2  = amin/NRESHI;

Nx2  = ceil(Sx/dx2);
Ny2  = ceil(Sy/dy2);

xa2  = linspace(xa(1), xa(Nx), Nx2);
ya2  = linspace(ya(1), ya(Ny), Ny2);

dx2  = xa2(2) - xa2(1);
dy2  = ya2(2) - ya2(1);

[Y2, X2] = meshgrid(ya2, xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD A GRAYSCALE UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL GRID
dxh = a/Nxh;
dyh = a/Nyh;
xah = [0:Nxh - 1]*dxh; xah = xah - mean(xah);
yah = [0:Nyh - 1]*dyh; yah = yah - mean(yah);
[YH, XH] = meshgrid(yah, xah);

% DEFINE VERTICES OF THE TRIANGLE
w  = 0.9*a;
h  = w*sqrt(3)/2;
v1 = [0 ; h/2];
v2 = [-w/2 ; -h/2];
v3 = [+w/2 ; -h/2];

% SECTION 1
p1 = v1;
p2 = v2;

D1 = (p2(2) - p1(2))*XH - (p2(1) - p1(1))*YH + p2(1)*p1(2) - p2(2)*p1(1);
D1 = -D1./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% SECTION 2
p1 = v1;
p2 = v3;

D2 = (p2(2) - p1(2))*XH - (p2(1) - p1(1))*YH + p2(1)*p1(2) - p2(2)*p1(1);
D2 = +D2./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% SECTION 3
p1 = v2;
p2 = v3;

D3 = (p2(2) - p1(2))*XH - (p2(1) - p1(1))*YH + p2(1)*p1(2) - p2(2)*p1(1);
D3 = -D3./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% BUILD THE UNIT CELL
UC = min(D1, D2);
UC = min(UC, D3);
UC = UC - min(UC(:));
UC = UC / max(UC(:));

if plotControl
    % SHOW UNIT CELL
    subplot(2,3,1);
    hh = imagesc(xah, yah, UC');
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('UNIT CELL');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FILL FRACTION SWEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE THRESHOLD DATA
gth_dat = linspace(0,1, 100);

% CALCULATE FILL FRACTION DATA
ff_dat = 0*gth_dat;

for n = 1 : length(gth_dat)
    % Generate Binary Unit Cell
    UCB = UC > gth_dat(n);
    
    % Calculate Fill Fraction
    ff_dat(n) = sum(UCB(:))/(Nxh*Nyh);  
end

if plotControl
    % SHOW FILL FRACTION SWEEP
    subplot(2, 3, 2);
    plot(gth_dat(1:n), ff_dat(1:n));
    xlim([gth_dat(1) max(gth_dat)]);
    ylim([0 1]);
    xlabel('Threshold Value');
    ylabel('Fill Fraction');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE FOURIER EXPANSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE FFT2
A = fftshift(fft2(UC))/(Nxh*Nyh);

% TRUNCATE FOURIER EXPANSION
p0  = ceil(Nxh/2);
q0  = ceil(Nyh/2);
np1 = p0 - floor(NP/2);
np2 = p0 + floor(NP/2);
nq1 = q0 - floor(NQ/2);
nq2 = q0 + floor(NQ/2);

AT = A(np1:np2, nq1:nq2);

% CALCULATE FOURIER AXES
pa = [-floor(NP/2):+floor(NP/2)];
qa = [-floor(NQ/2):+floor(NQ/2)];

% GRATING VECTOR EXPANSION
KX = 2*pi*pa/a;
KY = 2*pi*qa/a;
[KY, KX] = meshgrid(KY, KX);

if plotControl
    % SHOW TRUNCATED FFT2
    subplot(2,3,3);
    hh = imagesc(xah, yah, real(AT'));
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('TRUNCATED FFT2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERIOD (LOW RESOLUTION GRID)
PER = a*ones(Nx, Ny);

% ORIENTATION (LOW RESOLUTION GRID)
THETA = atan2(Y,X);

% FILL FACTOR (HIGH RESOLUTION GRID)
FF = (X2 - Sx/2).^2 + (Y2 - Sy/2).^2;
FF = FF - min(FF(:));
FF = FF / max(FF(:));
FF = 1 - FF;
FF = 0.05 + 0.3*FF;

if plotControl
    % SHOW INPUT DATA
    subplot(2,3,4);
    hh = imagesc(xah, yah, PER');
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('PER');

    subplot(2,3,5);
    hh = imagesc(xah, yah, THETA'/degrees);
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('THETA');

    subplot(2,3,6);
    hh = imagesc(xah, yah, FF');
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('FF');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE SPATIALLY-VARIATN LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATE LIST OF PLANAR GRATINGS
KLIST = [KX(:)' ; KY(:)'];
CLIST = AT(:);

% TRUNCATE LIST
cmax = max(abs(CLIST));
ind  = find(abs(CLIST) > cth*cmax);
CLIST = CLIST(ind);
KLIST = KLIST(:, ind);

%
% MAIN LOOP -- ITERATE OVER FOURIER EXPANSION
%

NK  = length(CLIST);
SVL = zeros(Nx2, Ny2);
for nk = 1 : NK
    
    % Calculate K-Function
    Kx         = KLIST(1, nk);
    Ky         = KLIST(2, nk);
    [theta, r] = cart2pol(Kx, Ky);
    [Kx, Ky]   = pol2cart(THETA+theta, r*a./PER);
    
    % Solve for Grating Phase
    PHI = svlsolve(Kx, Ky, dx, dy);
    
    % Interpolate to HIgh Resolution Grid
    PHI = interp2(ya, xa', PHI, ya2, xa2');
    
    % Calculate Analog Planar Grating
    GA = exp(1i*PHI);
    
    % Add Planar Grating to Overall Lattice
    SVL = SVL + CLIST(nk)*GA;
    
    % Refresh Graphics
    clf;
    
    if plotLatticeFormation
        % Show Planar Grating
        subplot(1,2,1);
        hh = imagesc(xah, yah, real(GA'));
        h2 = get(hh, 'Parent');
        set(h2, 'YDir', 'normal');
        axis equal tight;
        colorbar;
        title('PLANAR GRATING');

        % Show Planar Grating
        subplot(1,2,2);
        hh = imagesc(xah, yah, real(SVL'));
        h2 = get(hh, 'Parent');
        set(h2, 'YDir', 'normal');
        axis equal tight;
        colorbar;
        title('ANALOG LATTICE');

        % Force MATLAB to Draw Graphics
        if plot_now
            drawnow;
        end
        
        if MAKE_MOVIE == 1
            Frames = getframe(fig);
            writeVideo(vidObj, Frames);
        end
    end
    disp(['Displaying the Planar Grating ' num2str(nk) ' and the lattice.' ])
end

% CLEAN UP NUMERICAL NOISE
SVL = real(SVL);
SVL = SVL - min(SVL(:));
SVL = SVL / max(SVL(:));

% GENERATE BINARY LATTICE
GTH  = interp1(ff_dat, gth_dat, FF(:));
GTH  = reshape(GTH, Nx2, Ny2);
SVLB = SVL > GTH;

% SHOW FINAL LATTICE
clf;
hh = imagesc(xah, yah, SVLB');
h2 = get(hh, 'Parent');
set(h2, 'YDir', 'normal');
axis equal tight;
colorbar;
title('BINARY SPATIALLY-VARIANT LATTICE');


% SAVE BINARY LATTICE TO FILE
if save_lattice
    save svlattice2D SVLB dx2 dy2;
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






