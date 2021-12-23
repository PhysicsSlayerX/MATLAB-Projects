% rcwa3d_demo2.m
% Wavelength Sweep

% INITIALIZE MATLAB
close all;
clc; 
clear all;

% START TIMER
time1 = clock;

% UNITS
micrometers = 1;
nanometers  = 1e-3 * micrometers;
degrees     = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET FIGURE WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Color','k', 'Position', [0 42 1922 954]);
set(fig, 'Name', 'RCWA 3D');
set(fig, 'NumberTitle', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT CONTROLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_live = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
NLAM      = 200;
lam1      = 400 * nanometers;
lam2      = 800 * nanometers;
LAM0      = linspace(lam1, lam2, NLAM);
SRC.theta = 25*degrees;
SRC.phi   = 45*degrees;
SRC.pte   = 1;
SRC.ptm   = 0;

% DEVICE PARAMETERS
n_SiO = 1.4496;
n_SiN = 1.9360;
n_fs  = 1.5100;
a     = 1150/3 * nanometers;
r     =  500/3 * nanometers;
h1    =  230*2 * nanometers;
h2    =  345 * nanometers;
ax    =  0.7; 
ay    =  1;

DEV.er1 = 1.0;
DEV.ur1 = 1.0;
DEV.er2 = n_fs^2;
DEV.ur2 = 1.0;

DEV.t1 = [ a/2 ; -a*sqrt(3)/2 ];
DEV.t2 = [ a/2 ; +a*sqrt(3)/2 ];

% RCWA PARAMETERS
N1     = 512;
N2     = N1;
DEV.NP = 11;
DEV.NQ = DEV.NP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBLIQUE MESHGRID
p      = linspace(-0.5, +0.5, N1);
q      = linspace(-0.5, +0.5, N2);
[Q, P] = meshgrid(q, p);
XO     = P*DEV.t1(1) + Q*DEV.t2(1);
YO     = P*DEV.t1(2) + Q*DEV.t2(2);

% BUILD HEXAGONAL UNIT CELL
b   = a*sqrt(3);
RSQ = (XO/ax).^2 + ((YO - b/2)/ay).^2;
ER  = (RSQ <= r^2);
RSQ = (XO/ax).^2 + ((YO + b/2)/ay).^2;
ER  = ER | (RSQ <= r^2);
RSQ = ((XO - a/2)/ax).^2 + (YO/ay).^2;
ER  = ER | (RSQ <= r^2);
RSQ = ((XO + a/2)/ax).^2 + (YO/ay).^2;
ER  = ER | (RSQ <= r^2);

% CONVERT TO REAL MATERIALS
ER      = 1 - ER;
ERR     = 1 + (n_SiO^2 - 1)*ER;
URR     = ones(N1, N2);
DEV.L   = h1;

% ADD ADDITIONAL LAYERS
ERR(:,:,2) = (n_SiN^2)*ones(N1, N2);
URR(:,:,2) = ones(N1, N2);
DEV.L      = [DEV.L h2];

% SHOW DEVICE LAYERS
subplot(2,2,1);
pcolor(XO, YO, ERR(:,:,1));
shading interp;
axis equal tight;
colorbar;

subplot(2,2,2);
pcolor(XO, YO, ERR(:,:,2));
shading interp;
axis equal tight;
colorbar;


% COMPUTE CONVOLUTION MATRICES
NLAY = length(DEV.L);
NH   = DEV.NP*DEV.NQ;
DEV.ER   = zeros(NH,NH, NLAY);
DEV.UR   = zeros(NH,NH, NLAY);

for nlay   = 1 : NLAY
    DEV.ER(:,:,nlay) = convmat(ERR(:,:,nlay),DEV.NP,DEV.NQ);
    DEV.UR(:,:,nlay) = convmat(URR(:,:,nlay),DEV.NP,DEV.NQ);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALL RCWA3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE DATA RECORDS
REF = zeros(1, NLAM);
TRN = zeros(1, NLAM);

%
% MAIN LOOP -- WAVELENGTH SWEEP
% 
for nlam = 1 : NLAM
    
    % Call RCWA3D Function
    SRC.lam0 = LAM0(nlam);
    DAT      = rcwa3d(DEV, SRC);
    
    % Record Response
    REF(nlam) = DAT.REF;
    TRN(nlam) = DAT.TRN;
    CON = REF + TRN;
    
    if plot_live
        % Show Results
        subplot(2,2,[3 4]);
        plot(LAM0(1:nlam)/nanometers, CON(1:nlam), '--k','LineWidth', 2); hold on;
        plot(LAM0(1:nlam)/nanometers, REF(1:nlam), '-r', 'LineWidth', 2);
        plot(LAM0(1:nlam)/nanometers, TRN(1:nlam), '-b', 'LineWidth', 2); hold off;

        xlim([LAM0(1) LAM0(NLAM)]/nanometers);
        ylim([-0.05 1.05]);

        xlabel('Wavelength ($\lambda$)', 'Interpreter', 'LaTeX', 'FontSize', 16);
        ylabel('Responses' , 'Interpreter', 'LaTeX', 'FontSize', 16);

        title(['Optical Response at ' num2str(SRC.lam0/nanometers) ...
                ' nanometers'],'Interpreter', 'LaTeX', 'FontSize', 16);
        plot_darkmode;
        
        disp(['Computing at ' num2str(SRC.lam0/nanometers) ' nm']);
    end
    
    drawnow;
end

% PLOT THE ENTIRE RESPONSE
subplot(2,2,[3 4]);
plot(LAM0(1:NLAM)/nanometers, CON(1:NLAM), '--k','LineWidth', 2); hold on;
plot(LAM0(1:NLAM)/nanometers, REF(1:NLAM), '-r', 'LineWidth', 2);
plot(LAM0(1:NLAM)/nanometers, TRN(1:NLAM), '-b', 'LineWidth', 2); hold off;

xlim([LAM0(1) LAM0(NLAM)]/nanometers);
ylim([-0.05 1.05]);

xlabel('Wavelength ($\lambda$)', 'Interpreter', 'LaTeX', 'FontSize', 16);
ylabel('Responses' , 'Interpreter', 'LaTeX', 'FontSize', 16);

title(['Optical Responses'],'Interpreter', 'LaTeX', 'FontSize', 16);
plot_darkmode;

% REPORT RESULTS
disp(['Optical Response at ' num2str(LAM0(NLAM)/nanometers) ' nm' ] );
disp(['REF = ' num2str(100*REF(NLAM), '%6.2f') '%']);
disp(['TRN = ' num2str(100*TRN(NLAM), '%6.2f') '%']);
disp('================');
disp(['CON = ' num2str(100*CON(NLAM), '%6.2f') '%']);

% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);










