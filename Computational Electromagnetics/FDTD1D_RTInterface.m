
% BASED from the University of Texas - El Paso EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This code was written by Francis S. Dela Cruz, and was heavily referenced
% from the lectures on CEM (Computational Electromagnetics) by Dr. Raymond Rumpf. 
% The code was built from scratch with just the help of CEM Lectures.
% The basis of this code is from the Maxwell's Equation that was translated
% into matrices to be solved in MATLAB. This is a benchmarking code for
% future FDTD1D Calculations. Note that this was written with just a basic
% knowledge in MATLAB. Optimization and refactoring of codes are necessary
% to keep the runtime lower and simulation speed faster, as the coder
% implements better code syntax and alogrithms.

% Benchmarked as of August 5, 2020

% This MATLAB Program implements the Finite Difference Time Domain Method in 1D (FDTD-1D).

close all;
clc;
clear all;


t1 = clock;
opengl('software');
centimeters = 1;
millimeters = 0.1*centimeters;
meters = 100*centimeters;
seconds = 1;
minutes = 60*seconds;
c0 = 299792458*meters/seconds;
hertz = 1/seconds;
gigahertz = 1e9*hertz;

%% Simulation Parameters

freq = 1*gigahertz;
fmax = max(freq);
lambda_0 = c0/freq;
er_r = 1;
ur_r = 1;
er_src = 1;
ur_src = 1;
er_t = 1;
ur_t = 1;
ur_bc = 1;
er_bc = 1;
% Material Parameters
ur = 2.0;
er = 9.0;

L = 100;
h = 100;
th = 30.48*centimeters;

NRES = [10 4];
BUFZ = 0.5*[lambda_0 lambda_0];

numgraph_row = 2;
numgraph_col = 1;

freq_fine = 1; % make the frequency resolution finer
counter = 1; % Timestep controller - just for tests

% Movie Creation
MAKE_MOVIE = 1;
movie_title = 'FDTD1D_RTInterface.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    vidObj.FrameRate = 30;
    VidObj.Quality = 100;
    open(vidObj);
end

%% Consideration in Wavelength and Mechanical Resolution

% Wavelength
nref = sqrt(ur_r*er_r);
ntrn = sqrt(ur_t*er_t);
ngrid = sqrt(ur_bc*er_bc);
ndev = sqrt(ur*er);
nsrc = sqrt(ur_src*er_src);

nmax = max([nref ntrn ngrid ndev nsrc]);
lambda_min = c0/(fmax*nmax);
delta_lambda = lambda_min/NRES(1);

% Mechanical Features
mech_featX = [L];
mech_featY = [h];
mech_featZ = [th];

dminX = min(mech_featX);
dminY = min(mech_featY);
dminZ = min(mech_featZ);
dmin = min([dminX dminY dminZ]);
delta_d = dmin/NRES(2);

delta_x = min([delta_lambda delta_d]);
delta_y = delta_x;
delta_z = delta_x;

Mx = ceil(dminX/delta_x);
My = ceil(dminY/delta_y);
Mz = ceil(dminZ/delta_z);

dx = dminX/Mx;
dy = dminY/My;
dz = dminZ/Mz;

Nz = round(th/dz) + ceil((2*BUFZ(1) + BUFZ(2))/(nref*dz)) + 3; 

%% Initialize Materials to Free Space

nz1 = 2 + ceil(BUFZ(1)/dz) + 1;
nz2 = nz1 + round(th/dz) - 1;
nz3 = nz2 + 1;
nz4 = nz3 + round(th/2/dz) - 1;
nz5 = nz4 + 1;
nz6 = nz5 + round(th/5/dz) - 1;
ER = ones(1, Nz);
ER(nz1:nz2) = er;
ER(nz3:nz4) = er_r;
ER(nz5:nz6) = er*2;
UR = ones(1, Nz);
UR(nz1:nz2) = ur;
UR(nz3:nz4) = ur_r;
UR(nz5:nz6) = ur*2;
%% Compute Time Step
delta_min = min([dx dy dz]);
dt = (ngrid*delta_min)/(2*c0);


%% Compute for Update Coefficients

mEy = (c0*dt)./ER;
mHx = (c0*dt)./UR;

%% Field Initialization

Ey = zeros(1, Nz);
Hx = zeros(1, Nz);

%% Compute Source
nzref = 1;
nzsrc = nzref+1;
nztrn = Nz;

tau = 0.5/fmax;
time0 = 6*tau;

Tprop = (nmax*Nz*dz)/c0;
Ttotal = 12*tau + 5*Tprop;
STEPS = ceil(Ttotal/dt);

Time_element = [0:STEPS-1]*dt;
Eysrc  = exp(-((Time_element - time0)/tau).^2);
A = -sqrt(er_src/ur_src);
sigma_t = (nsrc*dz)/(2*c0)+dt/2;
Hxsrc = A.*exp(-((Time_element - time0 + sigma_t)/tau).^2);

% Initialize Fourier Transform
if freq_fine == 1
    NFREQ = STEPS;
else
    NFREQ = 100;
end
FREQ = linspace(0, freq, NFREQ);
K = exp(-1i*2*pi*dt.*FREQ);
REFe = zeros(1, NFREQ);
TRNe = zeros(1, NFREQ);
REFh = zeros(1, NFREQ);
TRNh = zeros(1, NFREQ);
SRC = zeros(1, NFREQ);

%% Initialize Boundary Terms

H2 = 0; H1 = H2; 
E2 = 0; E1 = E2;

%% Main FDTD Loop
fig1 = figure('Color','w');
set(fig1, 'Position', [1 41 1920 963]);
set(fig1, 'Name', 'FDTD-1D Analysis');
set(fig1, 'NumberTitle', 'off');

for T = 1:counter:STEPS
    % Update H from E (Perfect BC)
    H2 = H1; H1 = Hx(1);
    for nz = 1:Nz-1
        Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz + 1) - Ey(nz))/dz;
    end  
    Hx(Nz) = Hx(Nz) + mHx(Nz)*(E2 - Ey(Nz))/dz;
    % H-field Source
    Hx(nzsrc-1) = Hx(nzsrc-1) - (mHx(nzsrc-1)/dz)*Eysrc(T); 
    % Update E from H (Perfect BC)
    E2=E1; E1 = Ey(Nz);
    Ey(1) = Ey(1) + mEy(1)*(Hx(1) - H2)/dz;
    for nz = 2:Nz
        Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz - 1))/dz;
    end
    % E-Field Source
    Ey(nzsrc) = Ey(nzsrc) - (mEy(nzsrc)/dz)*Hxsrc(T);
    
    %Update Fourier Transforms
    for nf = 1: NFREQ
        REFe(nf) = REFe(nf) + (K(nf)^T)*Ey(nzref);
        TRNe(nf) = TRNe(nf) + (K(nf)^T)*Ey(nztrn);
        REFh(nf) = REFh(nf) + (K(nf)^T)*Hx(nzref);
        TRNh(nf) = TRNh(nf) + (K(nf)^T)*Hx(nztrn);
        SRC(nf) = SRC(nf) + (K(nf)^T)*Eysrc(T);
    end
    if MAKE_MOVIE == 1
        clf;
    end   
    % Update Visualization
    subplot(numgraph_row, numgraph_col, 1);
    x_axis = [0:Nz-1]*dz;
    x_fill1 = [nz1 nz2 nz2 nz1 nz1]*dz;
    x_fill2 = [nz4 nz6 nz6 nz4 nz4]*dz;
    y_fill = [-3 -3 3 3 -3];
    fill(x_fill1,y_fill,'g'); hold on;
    fill(x_fill2,y_fill,'g');
    axis([min(x_axis) max(x_axis) -3 3]);
    plot(x_axis, Ey,'-b', 'LineWidth', 3);
    plot(x_axis, Hx,'-r', 'LineWidth', 3); 
    hold off;
    axis([min(x_axis) max(x_axis) -3 3]);
    xlabel('z (cm)');
    ylabel('Field Amplitude','Rotation',90);
    title(['Timestep ' num2str(T) ' of ' num2str(ceil(Ttotal/dt)-1)]);
    h = legend('Device','Ey','Hx','Location','NorthEastOutside');
    set(h,'LineWidth',2);
    
    subplot(numgraph_row, numgraph_col, 2);
    x_axis = FREQ;
    plot(x_axis, (abs(REFe./SRC).^2),'-r','LineWidth', 3); hold on;
    plot(x_axis, (abs(TRNe./SRC).^2),'-b', 'LineWidth', 3);
    plot(x_axis, (abs(REFe./SRC).^2)+(abs(TRNe./SRC).^2),'--k', 'LineWidth', 3); hold off;
    axis([min(FREQ) max(FREQ) 0 1.25]);
    xlabel('Frequency (GHz)');
    ylabel('Reflectance and Transmittance','Rotation',90);
    title('Spectral Response - Ey');
    h = legend('REFe','TRNe','CONe','Location','NorthEastOutside');
    set(h,'LineWidth',2);
    drawnow;
    % Add Frame to AVI
    if MAKE_MOVIE == 1
        Frames = getframe(fig1);
        writeVideo(vidObj, Frames);
    end
end
if MAKE_MOVIE == 1
        close(vidObj);
end
% Finish Fourier Transform
REFe = REFe*dt;
TRNe = TRNe*dt;
REFh = REFh*dt;
TRNh = TRNh*dt;
SRC = SRC*dt;

% Compute Reflectance and Transmittance

REFe = abs(REFe./SRC).^2;
TRNe = abs(TRNe./SRC).^2;
CONe = REFe + TRNe;

REFh = abs(REFh./SRC).^2;
TRNh = abs(TRNh./SRC).^2;
CONh = REFh + TRNh;

RdBe = 10*log10(REFe);
TdBe = 10*log10(TRNe);
CdBe = 10*log10(CONe);

RdBh = 10*log10(REFh);
TdBh = 10*log10(TRNh);
CdBh = 10*log10(CONh);

% Visualization
fig2 = figure('Color','w');
set(fig2, 'Position', [1 41 1920 963]);
set(fig2, 'Name', 'FDTD-1D Analysis - Source Function');
set(fig2, 'NumberTitle', 'off');
x_axis = Time_element;
plot(x_axis, Eysrc,'-b','LineWidth', 3); hold on;
plot(x_axis, Hxsrc,'-r','LineWidth', 3); hold off;
axis([min(Time_element) max(Time_element) -1 1]);
xlabel('Time (s)');
ylabel('Amplitude','Rotation',90);
title('Source Functions');
h = legend('Ey','Hx','Location','NorthEastOutside');
set(h,'LineWidth',2);
 
fig3 = figure('Color','w');
set(fig3, 'Position', [1 41 1920 963]);
set(fig3, 'Name', 'FDTD-1D Analysis - Spectral Response');
set(fig3, 'NumberTitle', 'off');

subplot(2, 2, 1);

x_axis = FREQ;
plot(x_axis, REFe,'-r', 'LineWidth', 2.5); hold on;
plot(x_axis, TRNe,'-b','LineWidth', 2.5);
plot(x_axis, CONe,'--k','LineWidth', 2.5); hold off;
axis([min(FREQ) max(FREQ) 0 1.25]);
xlabel('Frequency (GHz)');
ylabel('Reflectance and Transmittance','Rotation',90);
title('Spectral Response - Ey');
h = legend('REFe','TRNe','CONe','Location','NorthEastOutside');
set(h,'LineWidth',2);

subplot(2, 2, 2);

x_axis = FREQ;
plot(x_axis, RdBe,'-r','LineWidth', 2.5); hold on;
plot(x_axis, TdBe,'-b','LineWidth', 2.5);
plot(x_axis, CdBe,'--k','LineWidth', 2.5); hold off;
axis([min(FREQ) max(FREQ) -30 0]);
xlabel('Frequency (GHz)');
ylabel('Reflectance and Transmittance (dB)','Rotation',90);
title('Spectral Response - Ey(dB)');
h = legend('REFe(dB)','TRNe(dB)','CONe(dB)','Location','NorthEastOutside');
set(h,'LineWidth',2);

subplot(2, 2, 3);

x_axis = FREQ;
plot(x_axis, REFh,'-r','LineWidth', 2.5); hold on;
plot(x_axis, TRNh,'-b','LineWidth', 2.5);
plot(x_axis, CONh,'--k','LineWidth', 2.5); hold off;
axis([min(FREQ) max(FREQ) 0 1.25]);
xlabel('Frequency (GHz)');
ylabel('Reflectance and Transmittance','Rotation',90);
title('Spectral Response - Hx');
h = legend('REFh','TRNh','CONh','Location','NorthEastOutside');
set(h,'LineWidth',2);

subplot(2, 2, 4);

x_axis = FREQ;
plot(x_axis, RdBh,'-r','LineWidth', 2.5); hold on;
plot(x_axis, TdBh,'-b','LineWidth', 2.5);
plot(x_axis, CdBh,'--k','LineWidth', 2.5); hold off;
axis([min(FREQ) max(FREQ) -30 0]);
xlabel('Frequency (GHz)');
ylabel('Reflectance and Transmittance (dB)','Rotation',90);
title('Spectral Response - Hx(dB)');
h = legend('REFh(dB)','TRNh(dB)','CONh(dB)','Location','NorthEastOutside');
set(h,'LineWidth',2);

t2 = clock;
t = etime(t2,t1);
disp(['Elapsed time is ' num2str(t/minutes) ...
    ' minutes.']);

