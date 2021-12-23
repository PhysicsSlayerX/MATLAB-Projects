% rcwa3d_demo.m

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
SRC.lam0  = 1540 * nanometers;
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
DEV.NP = 3;
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

% CALL FUNCTION
DAT = rcwa3d(DEV, SRC);

CON = DAT.REF + DAT.TRN;

% REPORT RESULTS
disp(['REF = ' num2str(100*DAT.REF, '%6.2f') '%']);
disp(['TRN = ' num2str(100*DAT.TRN, '%6.2f') '%']);
disp('================');
disp(['CON = ' num2str(100*CON, '%6.2f') '%']);









