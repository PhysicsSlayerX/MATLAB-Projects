% BASED from the University of Texas - El Paso EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This code was written by Francis S. Dela Cruz, and was heavily referenced
% from the lectures on CEM (Computational Electromagnetics). The code was
% built from scratch without any reference to prewritten codes available.
% The basis of this code is from the Maxwell's Equation that was translated
% into matrices to be solved in MATLAB. This is a benchmarking code for
% future TMM Calculations. Note that this was written with just a basic
% knowledge in MATLAB. Optimization and refactoring of codes are necessary
% to keep the runtime lower and simulation speed faster, as the coder
% implements better code algorithms and better code syntax. 

% Benchmarked as of April 24, 2020

% This MATLAB Program implements the Transfer Matrix Method (TMM).

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
lam0 = 2.7;                     %free space wavelength
theta = 57 * degrees;           %elevation angle
phi = 23 * degrees;             %azimuthal angle
pte = 1/sqrt(2);                %amplitude of TE polarization
ptm = 1i/sqrt(2);               %amplitude of TM polarization

% EXTERNAL MATERIALS
ur1 = 1.2;                      %permeability in the reflection region
er1 = 1.4;                      %permittivity in the reflection region
ur2 = 1.6;                      %permeability in the transmission region
er2 = 1.8;                      %permittivity in the transmission region

% DEFINE LAYERS
UR = [ 1 3 ];                   %array of permeabilities in each layer
ER = [ 2 1 ];                   %array of permittivities in each layer
L = [ 0.25 0.5 ];               %array of thickness of each layer

% DEFINE Identity Matrix and Zero Matrix

I = eye(2);
ZeroMat = zeros(2);

% Direction Normal to the Slab

n_hat = [0 0 -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the k0

k0 = 2*pi/lam0;

% Refractive Index in the reflection region and transmission region
n_inc = sqrt(ur1*er1);
n_trn = sqrt(ur2*er2);

% Calculate Longitudinal Component of the Incident Wave Vector 
kz_n = n_inc*cos(theta);

% Calculate Transverse Normalized Incident Wave Vectors
kx_n = n_inc*sin(theta)*cos(phi);
ky_n = n_inc*sin(theta)*sin(phi);

% Calculate Longitudinal Component of the Transmitted Wave Vector 
kz_nt = sqrt(ur2*er2 - kx_n^2 - ky_n^2);

% Calculate Longitudinal Component of the Wave Vector in the Gap Medium
urg = 1.0;                              %permeability in the gap medium
erg = 1.0 + kx_n^2 + ky_n^2;            %permittivity in the gap medium
kz_n_gap = sqrt(urg*erg - kx_n^2 - ky_n^2);

% Calculate the Gap Medium Parameters                        

Qg = [kx_n*ky_n 1+ky_n^2; -(1+kx_n^2) -kx_n*ky_n];
Vg = -1i*Qg;

% Initialize Global S-Matrix with Reflection Side S-Matrix (Add the
% reflection side S-matrix here)
SG11 = ZeroMat;
SG12 = I;
SG21 = I;
SG22 = ZeroMat;

S_global = [SG11 SG12; SG21 SG22];

% Calculate Parameter for ith Layer (revised this)
kz_i = zeros(1, length(L));
Qi = cell(1, length(L));               %preallocation - reducing the memory intake
OMEGA_i = cell(1, length(L));
Vi = cell(1, length(L));
for i = 1:length(L)
    Qi{i} = zeros(2);
    OMEGA_i{i} = zeros(2);
    Vi{i} = zeros(2);
end

for i = 1:length(L)
    kz_i(i) = sqrt(UR(i)*ER(i) - kx_n^2 - ky_n^2);
    Qi{i} = (1/UR(i))*[kx_n*ky_n, UR(i)*ER(i)-kx_n^2; ky_n^2-UR(i)*ER(i), -kx_n*ky_n];        
    OMEGA_i{i} = 1i*kz_i(i)*I;
    Vi{i} = Qi{i}/(OMEGA_i{i});
end

% Calculate S-Matrix for the ith Layer
Ai = cell(1, length(L));                  %preallocation
Bi = cell(1, length(L));
Xi = cell(1, length(L));
Di = cell(1, length(L));
S11_i = cell(1, length(L));
S12_i = cell(1, length(L));

for i = 1:length(L)
    Ai{i} = zeros(2);                     %Filling up the cells with 2x2 zero matrices
    Bi{i} = zeros(2);
    Xi{i} = zeros(2);
    Di{i} = zeros(2);
    S11_i{i} = zeros(2);
    S12_i{i} = zeros(2);
end
% Algorithm for getting the S-matrix components
for i = 1:length(L)
    Ai{i} = I + Vi{i}\Vg;
    Bi{i} = I - Vi{i}\Vg;
    Xi{i} = expm(OMEGA_i{i}*k0*L(i)*lam0);
    Di{i} = Ai{i} - Xi{i}*Bi{i}/Ai{i}*Xi{i}*Bi{i};
    S11_i{i} = Di{i}\(Xi{i}*Bi{i}/Ai{i}*Xi{i}*Ai{i} - Bi{i});
    S12_i{i} = Di{i}\Xi{i}*(Ai{i}-Bi{i}/Ai{i}*Bi{i});
end

S21_i = S12_i;
S22_i = S11_i;

% Update Global S-Matrix

D_i = cell(1, length(L));
Fi = cell(1, length(L));

for i = 1:length(L)
    D_i{i} = zeros(2);                     %Filling up the cells with 2x2 zero matrices
    Fi{i} = zeros(2);
end

D_i{1} = SG12/(I - S11_i{1}*SG22);
Fi{1} = S21_i{1}/(I - SG22*S11_i{1});

SG11 = SG11 + D_i{1}*S11_i{1}*SG21;
SG12 = D_i{1}*S12_i{1};
SG21 = Fi{1}*SG21;
SG22 = S22_i{1} + Fi{1}*SG22*S12_i{1};

D_i{2} = SG12/(I - S11_i{2}*SG22);
Fi{2} = S21_i{2}/(I - SG22*S11_i{2});

SG11 = SG11 + D_i{2}*S11_i{2}*SG21;
SG12 = D_i{2}*S12_i{2};
SG21 = Fi{2}*SG21;
SG22 = S22_i{2} + Fi{2}*SG22*S12_i{2};

% Connect to External Regions - REF and TRN Regions

% REF Region
kz_n_ref = sqrt(ur1*er1 - kx_n^2 - ky_n^2);
Qref = (1/ur1)*[kx_n*ky_n, ur1*er1-kx_n^2; ky_n^2-ur1*er1, -kx_n*ky_n];
OMEGAref = 1i*kz_n_ref*I;
Vref = Qref/OMEGAref;

Aref = I + Vg\Vref;
Bref = I - Vg\Vref;

SR11 = -Aref\Bref;
SR12 = 2*I/Aref;
SR21 = 0.5*I*(Aref - Bref/Aref*Bref);
SR22 = Bref/Aref;

% TRN Region

kz_n_trn = sqrt(ur2*er2 - kx_n^2 - ky_n^2);
Qtrn = (1/ur2)*[kx_n*ky_n, ur2*er2-kx_n^2; ky_n^2-ur2*er2, -kx_n*ky_n];
OMEGAtrn = 1i*kz_n_trn*I;
Vtrn = Qtrn/OMEGAtrn;

Atrn = I + Vg\Vtrn;
Btrn = I - Vg\Vtrn;

ST11 = Btrn/Atrn;
ST12 = 0.5*I*(Atrn - Btrn/Atrn*Btrn);
ST21 = 2*I/Atrn;
ST22 = -Atrn\Btrn;

% Finalizing the Global S-matrix by connecting the reflection and the
% transmission region to the S-Matrix of the device.

% Sglobal = Sglobal RSP S-trn

Dt = SG12/(I - ST11*SG22);
Ft = ST21/(I - SG22*ST11);

SG11 = SG11 + Dt*ST11*SG21;
SG12 = Dt*ST12;
SG21 = Ft*SG21;
SG22 = ST22 + Ft*SG22*ST12;

% Sglobal = S-ref RSP S-global

Dr = SR12/(I - SG11*SR22);
Fr = SG21/(I - SR22*SG11);

SG22 = SG22 + Fr*SR22*SG12;
SG21 = Fr*SR21;
SG12 = Dr*SG12;
SG11 = SR11 + Dr*SG11*SR21;

% Solving the Scatterring Problem

% Calculate the Source (Polarization components), then check for the |P| = 1

k_n_incident = [kx_n; ky_n; kz_n];

if k_n_incident(3,1) >= 1
    ate_hat = [0 1 0];
else
    ate_hat = cross(k_n_incident, n_hat)/norm(cross(k_n_incident, n_hat));
end

atm_hat = cross(ate_hat, k_n_incident)/norm(cross(ate_hat, k_n_incident));

P = pte*ate_hat + ptm*atm_hat;

normalized_P = norm(P); 

Esrcxy = [P(1);P(2)];

% Calculate for Transmitted and Reflected Fields

Erefxy = SG11*Esrcxy;

Etrnxy = SG21*Esrcxy;

% Calculate for the Longitudinal Field Components

Ezref = -(kx_n*Erefxy(1) + ky_n*Erefxy(2))/kz_n;
Eztrn = -(kx_n*Etrnxy(1) + ky_n*Etrnxy(2))/kz_nt;

% Calculate Transmittance and Reflectance

EREF = [Erefxy(1); Erefxy(2); Ezref];
ETRN = [Etrnxy(1); Etrnxy(2); Eztrn];

R = transpose(EREF)*conj(EREF);
T = (transpose(ETRN)*conj(ETRN))*(real(kz_nt/ur2)/real(kz_n/ur1));

% Check the conservation R+T = 1

CON = R + T;

% Displaying Texts and Values
disp(['R = ' num2str(R)]);
disp(['T = ' num2str(T)]);
disp(['CON = ' num2str(CON)]);

