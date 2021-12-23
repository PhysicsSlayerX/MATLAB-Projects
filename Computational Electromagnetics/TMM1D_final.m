% BASED from the University of Texas - El Paso EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This code was written by Francis S. Dela Cruz, and was heavily referenced
% from the lectures on CEM (Computational Electromagnetics). The code was
% built from scratch without any reference from prewritten codes available.
% The basis of this code is from the Maxwell's Equation that was translated
% into matrices to be solved in MATLAB. This is a benchmarking code for
% future TMM Calculations. Note that this was written with just a basic
% knowledge in MATLAB. Optimization and refactoring of codes are necessary
% to keep the runtime lower and simulation speed faster, as the coder
% implements better code syntax and alogrithms.

% Benchmarked as of May 11, 2020

% This MATLAB Program implements the Transfer Matrix Method (TMM).

% INITIALIZE MATLAB
close all;
clc;
clear all;

% START TIMER
t1 = clock;

% UNITS
degrees = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS

startPoint = 0.1;
stopPoint = 5.4;
numPoint = 10000;
LAMBDA = linspace(startPoint, stopPoint -1 , numPoint);

SRC.lam0 = 2.7;                         %free space wavelength
SRC.theta = 57 * degrees;               %elevation angle
SRC.phi = 23 * degrees;                 %azimuthal angle
SRC.pte = 1/sqrt(2);                    %amplitude of TE polarization
SRC.ptm = 1i/sqrt(2);                   %amplitude of TM polarization

% EXTERNAL MATERIALS
DEV.ur1 = 1.2;                          %permeability in the reflection region
DEV.er1 = 1.4;                          %permittivity in the reflection region
DEV.ur2 = 1.6;                          %permeability in the transmission region
DEV.er2 = 1.8;                          %permittivity in the transmission region

% DEFINE LAYERS
DEV.UR = [1 3];                        %array of permeabilities in each layer of unit cell
DEV.ER = [2 1];                        %array of permittivities in each layer of unit cell
DEV.L0 = [0.25 0.5];                    %array of thickness of each layer of the unit cell
DEV.lam0 = 2.7;
DEV.NP = 10;                            %number of times to cascade the unit cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocation of Output Parameters

NLAM = length(LAMBDA);
REF = zeros(1,NLAM);
TRN = zeros(1,NLAM);
CON = zeros(1, NLAM);

% CALCULATE DATA
for nlam = 1:NLAM
    SRC.lam0 = LAMBDA(nlam);
    DAT = tmm1d(DEV,SRC);
    REF(nlam) = DAT.REF;
    TRN(nlam) = DAT.TRN;
    CON(nlam) = DAT.CON;
end

% PLOTTING DATA
plotData(REF, TRN, CON, LAMBDA, NLAM);

% STOP TIMER
t2 = clock;
t = etime(t2,t1);
disp(['Elapsed time is ' num2str(t) ...
    ' seconds.']);

function DAT = tmm1d(DEV,SRC)
% TMM1D One-Dimensional Transfer Matrix Method
%
% DAT = tmm1d(DEV,SRC);
%
% Homework #5, Problem #3
% EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% INPUT ARGUMENTS
% ================
% DEV Device Parameters
%
% .er1 relative permittivity in reflection region
% .ur1 relative permeability in reflection region
% .er2 relative permittivity in transmission region
% .ur2 relative permeability in transmission region
%
% .ER array containing permittivity of each layer in unit cell
% .UR array containing permeability of each layer in unit cell
% .L array containing thickness of each layer in unit cell
% .NP number of unit cells to cascade
%
% SRC Source Parameters
%
% .lam0 free space wavelength
%
% .theta elevation angle of incidence (radians)
% .phi azimuthal angle of incidence (radians)
%
% .ate amplitude of TE polarization
% .atm amplitude of TM polarization
%
% OUTPUT ARGUMENTS
% ================
% DAT Output Data
%
% .REF Reflectance
% .TRN Transmittance
% Calculation of the k0

    k0 = 2*pi/SRC.lam0;
    DEV.L = DEV.L0*DEV.lam0;
% Refractive Index in the reflection region and transmission region
    n_inc = sqrt(DEV.ur1*DEV.er1);

% Calculate Transverse Normalized Incident Wave Vectors
    
    kx_n = n_inc*sin(SRC.theta)*cos(SRC.phi);
    ky_n = n_inc*sin(SRC.theta)*sin(SRC.phi);
    
% Calculate Longitudinal Component of the Incident and Transmitted Wave Vector 
    kz_n = sqrt(DEV.ur1*DEV.er1 - kx_n^2 - ky_n^2);
    kz_nt = sqrt(DEV.ur2*DEV.er2 - kx_n^2 - ky_n^2);

% Calculate Longitudinal Component of the Wave Vector in the Gap Medium
    %{
    urg = 1.0;                              %permeability in the gap medium
    erg = 1.0 + kx_n^2 + ky_n^2;            %permittivity in the gap medium
    kz_n_gap = sqrt(urg*erg - kx_n^2 - ky_n^2);
    %}

% Calculate the Gap Medium Parameters                        

    Qg = [kx_n*ky_n 1+ky_n^2; -(1+kx_n^2) -kx_n*ky_n];
    Vg = -1i*Qg;

% Initialize Global S-Matrix 
    I = eye(2);
    %{
    ZeroMat = zeros(2);
    SG.S11 = ZeroMat;
    SG.S12 = I;
    SG.S21 = I;
    SG.S22 = ZeroMat;
    %}
    
% Calculate Parameter for ith Layer (revised this)
    kz_i = zeros(1, length(DEV.L));
    Qi = cell(1, length(DEV.L));               %preallocation - reducing the memory intake
    OMEGA_i = cell(1, length(DEV.L));
    Vi = cell(1, length(DEV.L));
    for i = 1:length(DEV.L)
        Qi{i} = zeros(2);
        OMEGA_i{i} = zeros(2);
        Vi{i} = zeros(2);
    end

    for i = 1:length(DEV.L)
        kz_i(i) = sqrt(DEV.UR(i)*DEV.ER(i) - kx_n^2 - ky_n^2);
        Qi{i} = (1/DEV.UR(i))*[kx_n*ky_n, DEV.UR(i)*DEV.ER(i)-kx_n^2; ky_n^2-DEV.UR(i)*DEV.ER(i), -kx_n*ky_n];        
        OMEGA_i{i} = 1i*kz_i(i)*I;
        Vi{i} = Qi{i}/(OMEGA_i{i});
    end

% Calculate S-Matrix for the ith Layer
    Ai = cell(1, length(DEV.L));                  %preallocation
    Bi = cell(1, length(DEV.L));
    Xi = cell(1, length(DEV.L));
    Di = cell(1, length(DEV.L));
    S11_i = cell(1, length(DEV.L));
    S12_i = cell(1, length(DEV.L));

    for i = 1:length(DEV.L)
        Ai{i} = zeros(2);                     %Filling up the cells with 2x2 zero matrices
        Bi{i} = zeros(2);
        Xi{i} = zeros(2);
        Di{i} = zeros(2);
        S11_i{i} = zeros(2);
        S12_i{i} = zeros(2);
    end
% Algorithm for getting the S-matrix components
    for i = 1:length(DEV.L)
        Ai{i} = I + Vi{i}\Vg;
        Bi{i} = I - Vi{i}\Vg;
        Xi{i} = expm(OMEGA_i{i}*k0*DEV.L(i));
        Di{i} = Ai{i} - Xi{i}*Bi{i}/Ai{i}*Xi{i}*Bi{i};
        S11_i{i} = Di{i}\(Xi{i}*Bi{i}/Ai{i}*Xi{i}*Ai{i} - Bi{i});
        S12_i{i} = Di{i}\Xi{i}*(Ai{i}-Bi{i}/Ai{i}*Bi{i});
    end

    S21_i = S12_i;
    S22_i = S11_i;

% One Unit Cell S-Matrix

    if isempty(DEV.L)
        SD.S11 = zeros(2);
        SD.S12 = eye(2);
        SD.S21 = eye(2);
        SD.S22 = zeros(2);
    else
        SD.S11 = S11_i{1};
        SD.S12 = S12_i{1};
        SD.S21 = S21_i{1};
        SD.S22 = S22_i{1};

        for i = 2:length(DEV.L)
            SDn.S11 = S11_i{i};
            SDn.S12 = S12_i{i};
            SDn.S21 = S21_i{i};
            SDn.S22 = S22_i{i};
            SD = star(SD, SDn);
        end

    end

% Cascading and Doubling Algorithm

    SG = cascn(SD, DEV.NP);

% Connect to External Regions - REF and TRN Regions

    % REF Region
    kz_n_ref = sqrt(DEV.ur1*DEV.er1 - kx_n^2 - ky_n^2);
    Qref = (1/DEV.ur1)*[kx_n*ky_n, DEV.ur1*DEV.er1-kx_n^2; ky_n^2-DEV.ur1*DEV.er1, -kx_n*ky_n];
    OMEGAref = 1i*kz_n_ref*I;
    Vref = Qref/OMEGAref;

    Aref = I + Vg\Vref;
    Bref = I - Vg\Vref;

    SR11 = -Aref\Bref;
    SR12 = 2*I/Aref;
    SR21 = 0.5*I*(Aref - Bref/Aref*Bref);
    SR22 = Bref/Aref;

    % TRN Region

    kz_n_trn = sqrt(DEV.ur2*DEV.er2 - kx_n^2 - ky_n^2);
    Qtrn = (1/DEV.ur2)*[kx_n*ky_n, DEV.ur2*DEV.er2-kx_n^2; ky_n^2-DEV.ur2*DEV.er2, -kx_n*ky_n];
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

    Dt = SG.S12/(I - ST11*SG.S22);
    Ft = ST21/(I - SG.S22*ST11);

    SG.S11 = SG.S11 + Dt*ST11*SG.S21;
    SG.S12 = Dt*ST12;
    SG.S21 = Ft*SG.S21;
    SG.S22 = ST22 + Ft*SG.S22*ST12;

    % Sglobal = S-ref RSP S-global

    Dr = SR12/(I - SG.S11*SR22);
    Fr = SG.S21/(I - SR22*SG.S11);

    SG.S22 = SG.S22 + Fr*SR22*SG.S12;
    SG.S21 = Fr*SR21;
    SG.S12 = Dr*SG.S12;
    SG.S11 = SR11 + Dr*SG.S11*SR21;

% Solving the Scatterring Problem

% Calculate the Source (Polarization components), then check for the |P| = 1

    k_n_incident = k0*[kx_n ky_n kz_n];
    n_hat = [0 0 -1];

    if SRC.theta == 0
        ate_hat = [0 1 0];
    else
        ate_hat = cross(k_n_incident, n_hat)/norm(cross(k_n_incident, n_hat));
    end

    atm_hat = cross(ate_hat, k_n_incident)/norm(cross(ate_hat, k_n_incident));

    P = SRC.pte*ate_hat + SRC.ptm*atm_hat;

    normalized_P = norm(P); 

    Esrcxy = [P(1);P(2)];

    % Calculate for Transmitted and Reflected Fields

    Erefxy = SG.S11*Esrcxy;

    Etrnxy = SG.S21*Esrcxy;

    % Calculate for the Longitudinal Field Components

    Ezref = -(kx_n*Erefxy(1) + ky_n*Erefxy(2))/kz_n;
    Eztrn = -(kx_n*Etrnxy(1) + ky_n*Etrnxy(2))/kz_nt;

    % Calculate Transmittance and Reflectance

    EREF = [Erefxy(1); Erefxy(2); Ezref];
    ETRN = [Etrnxy(1); Etrnxy(2); Eztrn];

    DAT.REF = transpose(EREF)*conj(EREF);
    DAT.TRN = (transpose(ETRN)*conj(ETRN))*(real(kz_nt/DEV.ur2)/real(kz_n/DEV.ur1));

    % Check the conservation R+T = 1

    DAT.CON = DAT.REF + DAT.TRN;

    % Displaying Texts and Values
    disp(['R = ' num2str(DAT.REF)]);
    disp(['T = ' num2str(DAT.TRN)]);
    disp(['CON = ' num2str(DAT.CON)]);
    
end

function S = star(SA,SB)
    % STAR Redheffer Star Product
    %
    % S = star(SA,SB)
    %
    % INPUT ARGUMENTS
    % ================
    % SA First Scattering Matrix
    % .S11 S11 scattering parameter
    % .S12 S12 scattering parameter
    % .S21 S21 scattering parameter
    % .S22 S22 scattering parameter
    %
    % SB Second Scattering Matrix
    % .S11 S11 scattering parameter
    % .S12 S12 scattering parameter
    % .S21 S21 scattering parameter
    % .S22 S22 scattering parameter
    %
    % OUTPUT ARGUMENTS
    % ================
    % S Combined Scattering Matrix
    % .S11 S11 scattering parameter
    % .S12 S12 scattering parameter
    % .S21 S21 scattering parameter
    % .S22 S22 scattering parameter
    I = eye(size(SA.S11));
    D = SA.S12/(I - SB.S11*SA.S22);
    F = SB.S21/(I - SA.S22*SB.S11);
    S.S11 = SA.S11 + D*(SB.S11*SA.S21);
    S.S12 = D*SB.S12;
    S.S21 = F*SA.S21;
    S.S22 = SB.S22 + F*(SA.S22*SB.S12);
end

function SC = cascn(S, N)
    I = eye(size(S.S11));
    ZeroMat = zeros(size(S.S11));
    SC.S11 = ZeroMat;
    SC.S12 = I;
    SC.S21 = I;
    SC.S22 = ZeroMat;

    N_bin = de2bi(N);

    Sbin = S;

    for i = 1:length(N_bin)
        if N_bin(i)== 1
            SC = star(SC, Sbin);
            if i <= length(N_bin) - 1
                disp(num2str(i));
                Sbin = star(Sbin,Sbin);
            else
                break
            end
        elseif N_bin(i)== 0
            Sbin = star(Sbin,Sbin);
        end
    end
end

function [f1,f2,f3] = plotData(ref, trn, con, xAxis, numEL)
    f1 = ref;
    f2 = trn;
    f3 = con;
    lineWidth = 1.5;
    % OPEN FIGURE WINDOW
    F1 = figure('Color','w');
    % PLOT DATA
    h = plot(xAxis,f1,'-b','LineWidth',lineWidth);
    hold on;
    plot(xAxis,f2,'--r','LineWidth',lineWidth);
    hold on;
    plot(xAxis,f3,'-.k','LineWidth',lineWidth);
    hold off;
    % SET AXIS LIMITS
    xlim([0.5 xAxis(numEL)]);
    ylim([0 1.1]);
    % MAKE LINES THICK AND FONTS BIGGER
    h2 = get(h,'Parent');
    set(h2,'LineWidth',3,'FontSize',18);
    % SET TICK MARKS

    xm = [0:xAxis(numEL)];
    xt = {};
    for m = 1 : length(xm)
    xt{m} = num2str(xm(m),'%3.2f');
    end
    set(h2,'XTick',xm,'XTickLabel',xt);
    ym = [0:0.1+1];
    yt = {};
    for m = 1 : length(ym)
    yt{m} = num2str(ym(m),'%2.1f');
    end
    set(h2,'YTick',ym,'YTickLabel',yt);

    % LABEL AXES
    xlabel('$ \textrm{Wavelength, } \lambda $','Interpreter','latex');
    ylabel('$ \textrm{Response} $','Interpreter','latex');
    % ADD LEGEND
    h = legend('Wavelength(\lambda)','Response','Location','NorthOutside');
    set(h,'LineWidth',2);
end
