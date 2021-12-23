% rcwa3d_demo2.m
% Wavelength Sweep

% INITIALIZE MATLAB
close all;
clc; 
clear all;

% UNITS
micrometers = 1;
nanometers  = 1e-3 * micrometers;
degrees     = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
NLAM      = 200;
lam1      = 1530 * nanometers;
lam2      = 1550 * nanometers;
LAM0      = linspace(lam1, lam2, NLAM);
SRC.theta = 0*degrees;
SRC.phi   = 0*degrees;
SRC.pte   = 1;
SRC.ptm   = 0;

% DEVICE PARAMETERS
n_SiO = 1.4496;
n_SiN = 1.9360;
n_fs  = 1.5100;
a     = 1150 * nanometers;
r     =  400 * nanometers;
h1    =  230 * nanometers;
h2    =  345 * nanometers;

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
RSQ = XO.^2 + (YO - b/2).^2;
ER  = (RSQ <= r^2);
RSQ = XO.^2 + (YO + b/2).^2;
ER  = ER | (RSQ <= r^2);
RSQ = (XO - a/2).^2 + (YO).^2;
ER  = ER | (RSQ <= r^2);
RSQ = (XO + a/2).^2 + (YO).^2;
ER  = ER | (RSQ <= r^2);

% CONVERT TO REAL MATERIALS
ER  = 1 - ER;
ERR = 1 + (n_SiO^2 - 1)*ER;
URR = ones(N1, N2);
DEV.L   = h1;

% ADD ADDITIONAL LAYERS
ERR(:,:,2) = (n_SiN^2)*ones(N1, N2);
URR(:,:,2) = ones(N1, N2);
DEV.L          = [DEV.L h2];

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
    
    % Show Results
    plot(LAM0(1:nlam)/nanometers, CON(1:nlam), '--k'); hold on;
    plot(LAM0(1:nlam)/nanometers, REF(1:nlam), '-r');
    plot(LAM0(1:nlam)/nanometers, TRN(1:nlam), '-b'); hold off;
    
    xlim([LAM0(1) LAM0(NLAM)]/nanometers);
    ylim([-0.05 1.05]);
    
    drawnow;
end











