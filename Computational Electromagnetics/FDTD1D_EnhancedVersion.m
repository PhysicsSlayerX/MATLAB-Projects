
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
micrometers = 1;
nanometers = (1/1000)*micrometers;
meters = 1000000*micrometers;
seconds = 1;
minutes = 60*seconds;
c0 = 299792458*meters/seconds;
hertz = 1/seconds;
gigahertz = 1e9*hertz;
terahertz = 1e12*hertz;

%% Simulation Parameters

% freq = 5*gigahertz;
% fmax = max(freq);
% lambda_0 = c0/freq;
lambda_0 = 980*nanometers;
fmax = (c0/lambda_0);
freq = fmax;
er_r = 1;
ur_r = 1;
er_src = 1;
ur_src = 1;
er_t = 1;
ur_t = 1;
ur_bc = 1;
er_bc = 1;
% Material Parameters
ur = 1;
er = 1.5^2;
ur_rl = 1;
er_rl = 2.0^2;
ur_r2 = 1;
er_r2 = 5.0^2;

GratingPeriod = 5;

L = 100;
h = 100;
d1 = ((lambda_0)/(4*sqrt(ur*er)));
d2 = ((lambda_0)/(4*sqrt(ur_rl*er_rl)));
d3 = ((lambda_0)/(4*sqrt(ur_r2*er_r2)));

% For time step scaling

C1 = [20 10];

NRES = [10 4];
BUFZ = 0.50*[lambda_0 lambda_0];

numgraph_row = 2;
numgraph_col = 1;

freq_fine = 1; % make the frequency resolution finer
counter = 1; % Timestep controller - just for tests

gaussiansource = 1;
sinesource = 0;

% Movie Creation
MAKE_MOVIE = 1;

movie_title = 'FDTD1D_EV2.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    vidObj.FrameRate = 30;
    VidObj.Quality = 100;
    open(vidObj);
end

%% Consideration in Wavelength and Mechanical Resolution

% Wavelength
nref    = sqrt(ur_r*er_r);
ntrn    = sqrt(ur_t*er_t);
ngrid   = sqrt(ur_bc*er_bc);
nd1     = sqrt(ur*er);
nd2     = sqrt(ur_rl*er_rl);
nd3     = sqrt(ur_r2*er_r2);
nsrc    = sqrt(ur_src*er_src);

nmax = max([nref ntrn ngrid nd1 nsrc nd2]);
%lambda_min = c0/(fmax*nmax);
lambda_min = c0/(fmax*nmax);
delta_lambda = lambda_min/NRES(1);

% Mechanical Features
mech_featX = [L];
mech_featY = [h];
mech_featZ = [d1 d2 d3];

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

Nz = round((sum(mech_featZ)*GratingPeriod)/dz) + ceil((BUFZ(1) + BUFZ(2))/(nref*dz)) + 3; 

%% Initialize Materials to Free Space


ER = ones(1, Nz);

nz1 = 2 + ceil(0.5*BUFZ(1)/dz) + 1; %location of the top surface of the grating from the top of the grid
for build = 1:GratingPeriod
    b1 = nz1 + round(d1/dz) - 1; %location of the tooth base from the top of the grid
    c1 = b1+1; %location of the top surface of the grating substrate from the top of the grid
    dwow = c1 + round(d2/dz) - 1; %location of the bottom surface of the grating from the top of the grid
    e1 = dwow + 1;
    f1 = e1 + round(d3/dz) - 1;
    ER(nz1:b1) = er; %Generate and fill the first tooth with the appropriate dielectric values
    ER(c1:dwow) = er_rl; %Generate and fill the grating substrate
    ER(e1:f1) = er_r2; %Generate and fill the grating substrate
    nz1 = f1;
end

UR = ones(1, Nz);
for build = 1:GratingPeriod
    b1 = nz1 + round(d1/dz) - 1; %location of the tooth base from the top of the grid
    c1 = b1+1; %location of the top surface of the grating substrate from the top of the grid
    dwow = c1 + round(d2/dz) - 1; %location of the bottom surface of the grating from the top of the grid
    e1 = dwow + 1;
    f1 = e1 + round(d3/dz) - 1;
    UR(nz1:b1) = ur; %Generate and fill the first tooth with the appropriate dielectric values
    UR(c1:dwow) = ur_rl; %Generate and fill the grating substrate
    UR(e1:f1) = ur_r2; %Generate and fill the grating substrate
    nz1 = f1;
end

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

if gaussiansource == 1
    tau = 0.5/fmax;
    time0 = 6*tau;

    Tprop = (nmax*Nz*dz)/c0;
    Ttotal = C1(1)*tau + C1(2)*Tprop;
    STEPS = ceil(Ttotal/dt);

    Time_element = [0:STEPS-1]*dt;
    Eysrc  = exp(-((Time_element - time0)/tau).^2);
    A = -sqrt(er_src/ur_src);
    sigma_t = (nsrc*dz)/(2*c0)+dt/2;
    Hxsrc = A.*exp(-((Time_element - time0 + sigma_t)/tau).^2);
end

if sinesource == 1
    
    tau = 3/fmax;
    time0 = 3*tau;

    Tprop = (nmax*Nz*dz)/c0;
    Ttotal = C1(1)*tau + C1(2)*Tprop;
    STEPS = ceil(Ttotal/dt);

    Time_element = [0:STEPS-1]*dt;
    
    A_h = -sqrt(er_src/ur_src);
    sigma_t = (nsrc*dz)/(2*c0)+dt/2;
    Hxsrc(Time_element >= time0) = A_h.* sin(2*pi*fmax*Time_element(Time_element >= time0));
    Eysrc(Time_element >= time0) = sin(2*pi*fmax*Time_element(Time_element >= time0));
    Eysrc(Time_element < time0) = exp(-((Time_element(Time_element < time0) - time0)/tau).^2).* sin(2*pi*fmax*Time_element(Time_element < time0));
    Hxsrc(Time_element < time0) = A_h.*exp(-((Time_element(Time_element < time0) - time0 + sigma_t)/tau).^2).*sin(2*pi*fmax*Time_element(Time_element < time0));
    
end


% Initialize Fourier Transform
if freq_fine == 1
    NFREQ = STEPS;
else
    NFREQ = 100;
end
FREQ = linspace(0*freq, (freq*1.5), NFREQ);
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
    y_fill = [-3 -3 3 3 -3];
    
    nz1 = 2 + ceil(BUFZ(1)/dz) + 1;
    for build = 1:GratingPeriod
        b1 = nz1 + round(d1/dz) - 1; %location of the tooth base from the top of the grid
        c1 = b1+1; %location of the top surface of the grating substrate from the top of the grid
        dwow = c1 + round(d2/dz) - 1; %location of the bottom surface of the grating from the top of the grid
        e1 = dwow + 1;
        f1 = e1 + round(d3/dz) - 1;
        x_fill1 = [nz1 b1 b1 nz1 nz1]*dz;
        x_fill2 = [b1 dwow dwow b1 b1]*dz;
        x_fill3 = [dwow f1 f1 dwow dwow]*dz;
        fill(x_fill1,y_fill,'c'); hold on;
        fill(x_fill2,y_fill,'g');
        fill(x_fill3,y_fill,'m');
        nz1 = f1;
    end
    
    axis([min(x_axis) max(x_axis) -3 3]);
    plot(x_axis, Ey,'-b', 'LineWidth', 3);
    plot(x_axis, Hx,'-r', 'LineWidth', 3); 
    hold off;
    axis([min(x_axis) max(x_axis) -3 3]);
    xlabel('z (\mum)');
    ylabel('Field Amplitude','Rotation',90);
    title(['Timestep ' num2str(T) ' of ' num2str(ceil(Ttotal/dt)-1)]);
    h = legend('d1', 'd2','Ey','Hx','Location','NorthEastOutside');
    set(h,'LineWidth',2);
    
    subplot(numgraph_row, numgraph_col, 2);
    x_axis = FREQ;
    plot(x_axis, 10*log10(abs(REFe./SRC).^2),'-r','LineWidth', 3); hold on;
    plot(x_axis, 10*log10(abs(TRNe./SRC).^2),'-b', 'LineWidth', 3);
    plot(x_axis, 10*log10((abs(REFe./SRC).^2)+(abs(TRNe./SRC).^2)),'--k', 'LineWidth', 3); hold off;
    axis([min(x_axis) max(x_axis) -50 0.5]);
    xlabel('Frequency (Hz)');
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
plot(x_axis, Eysrc,'-b','LineWidth', 0.75); hold on;
plot(x_axis, Hxsrc,'-r','LineWidth', 0.75); hold off;
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

