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
set(fig, 'Name', 'Demonstration of 2D Fourier Expansion');
set(fig, 'NumberTitle', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MOVIE CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Movie Creation
MAKE_MOVIE = 1;
movie_title = '2D Fourier Expansion.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    vidObj.FrameRate = 30;
    VidObj.Quality = 100;
    open(vidObj);
end

% PLOT CONTROL
plotControl = 1;
plot_now = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL PARAMETERS
a      = 1;
w      = 0.9*a;

% GRID PARAMETERS
Nx = 255;
Ny = Nx;

% FOURIER EXPANSION PARAMETERS
NP = 21;
NQ = 21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID
dx = a/Nx;
dy = a/Ny;

xa = [0: Nx-1]*dx; xa = xa - mean(xa);
ya = [0: Ny-1]*dy; ya = ya - mean(ya);
[Y, X] = meshgrid(ya, xa);

% GRID INDICES OF TRIANGLE
h   = w*sqrt(3)/2;
ny  = round(h/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;

% BUILD TRIANGLE
UC = zeros(Nx, Ny);
for ny = ny1 : ny2
    f = 1 - (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f*w/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    UC(nx1:nx2, ny) = 1;
end

% SHOW UNIT CELL
subplot(2, 3, 1);
hh = imagesc(xa, ya, UC');
h2 = get(hh, 'Parent');
set(h2, 'YDir', 'normal');
axis equal tight;
colorbar;
title('UNIT CELL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE FOURIER EXPANSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE FFT2
A = fftshift(fft2(UC))/(Nx*Ny);

% TRUNCATE FFT2
p0 = ceil(Nx/2);
q0 = ceil(Ny/2);
np1 = p0 - floor(NP/2);
np2 = p0 + floor(NP/2);
nq1 = q0 - floor(NQ/2);
nq2 = q0 + floor(NQ/2);

AT = A(np1:np2, nq1:nq2);

% AXES FOR AT
pa = [-floor(NP/2): floor(NP/2)];
qa = [-floor(NQ/2): floor(NQ/2)];

% CALCULATE GRATING VECTOR EXPANSION
KX = 2*pi*pa/a;
KY = 2*pi*qa/a;
[KY, KX] = meshgrid(KY, KX);

if plotControl
    % SHOW FFT2
    subplot(2, 3, 2); 
    hh = imagesc(xa, ya, real(A'));
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('FFT2');

    % SHOW TRUNCATED FFT2
    subplot(2, 3, 3); 
    hh = imagesc(pa, qa, real(AT'));
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('TRUNCATED FFT2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECONSTRUCT UNIT CELL FROM FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUICK WAY
 A = zeros(Nx, Ny);
 A(np1:np2, nq1:nq2) = AT;
 UCT = real(ifft2(ifftshift(A)))*(Nx*Ny);

if plotControl
    % SHOW TRUNCATED FFT2
    subplot(2, 3, 3); 
    hh = imagesc(xa, ya, real(UCT'));
    h2 = get(hh, 'Parent');
    set(h2, 'YDir', 'normal');
    axis equal tight;
    colorbar;
    title('TRUNCATED UNIT CELL');
end

% SLOW WAY
UCT = zeros(Nx, Ny);
for nq = 1 : NQ
    for np = 1 : NP
        
        % Calculate Planar Grating
        G = exp(1i*(KX(np, nq)*X + KY(np, nq)*Y));
        
        % Add Grating to Overall Unit Cell
        UCT = UCT + AT(np, nq)*G;
        
        if plotControl 
        % Show Truncated FFT
            subplot(2,3,4);
            hh = imagesc(pa, qa, real(AT'));
            h2 = get(hh, 'Parent');
            set(h2, 'YDir', 'normal');
            axis equal tight;
            colorbar;
            title('TRUNCATED FFT2');
            hold on;
            x = pa(np) - 0.5 + [0 1 1 0 0];
            y = qa(nq) - 0.5 + [0 0 1 1 0];
            line(x, y, 'Color', 'k', 'LineWidth', 2);
            hold off;

            % Show Planar Grating
            subplot(2, 3, 5)
            hh = imagesc(xa, ya, real(G'));
            h2 = get(hh, 'Parent');
            set(h2, 'YDir', 'normal');
            axis equal tight;
            colorbar;
            title(['P = ' num2str(pa(np)) ', Q = ' ...
                    num2str(qa(nq) )])   

            % Show Reconstructed Unit Cell
            subplot(2, 3, 6)
            hh = imagesc(xa, ya, ifftshift(real(UCT')));
            h2 = get(hh, 'Parent');
            set(h2, 'YDir', 'normal');
            axis equal tight;
            colorbar;
            title('RECONSTRUCTED UNIT CELL') 

            % Force MATLAB to Draw Graphics 
            if plot_now
                drawnow;
            end
        end
        if MAKE_MOVIE == 1
            Frames = getframe(fig);
            writeVideo(vidObj, Frames);
        end
        disp(['Processing Planar Gratings: ' 'P = ' num2str(pa(np)) ', Q = ' ...
                    num2str(qa(nq) )])
    end
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








