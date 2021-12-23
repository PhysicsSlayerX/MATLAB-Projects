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
%Set Variable Precision
digits(4);
% START TIMER
t1 = clock;

% UNITS
degrees = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS


startPoint = 1;
stopPoint = 2;
numPoint = 100;
LAMBDA = linspace(startPoint, stopPoint , numPoint);

seed = 14;
rng(seed);                                  %Set the random seed (used seeds - 1 = train, 2 = test, 3 = sample)

SRC.lam0 = 1.0;                        %free space wavelength
SRC.theta = 0 * degrees;                %elevation angle
SRC.phi = 0 * degrees;                  %azimuthal angle
SRC.pte = 1;                            %amplitude of TE polarization
SRC.ptm = 0;                            %amplitude of TM polarization

% EXTERNAL MATERIALS
DEV.ur1 = 1.0;                          %permeability in the reflection region
DEV.er1 = 1.0;                          %permittivity in the reflection region
DEV.ur2 = 1.0;                          %permeability in the transmission region
DEV.er2 = 1.0;                          %permittivity in the transmission region

% DEFINE LAYERS
DEV.UR = [1 1 1 1];                         %array of permeabilities in each layer of unit cell
DEV_ER = (abs(1+(5-1).*rand(10000,4))).^2; %array of permittivities in each layer of unit cell
DEVThick = abs((1+(3-1).*rand(10000,4)));                   %array of thickness of each layer of the unit cell
DEV.lam0 = 1;
DEV.NP = 1;                            %number of times to cascade the unit cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocation of Output Parameters

NLAM = length(LAMBDA);
layerLength = size(DEVThick, 1);
REF = zeros(layerLength,NLAM);
TRN = zeros(layerLength,NLAM);
CON = zeros(layerLength, NLAM);

REFdB = zeros(layerLength,NLAM);
TRNdB = zeros(layerLength,NLAM);
CONdB = zeros(layerLength,NLAM);

% CALCULATE DATA
for layer = 1:layerLength
    for nlam = 1:NLAM
        SRC.lam0 = LAMBDA(nlam);
        DEV.L0 = DEVThick(layer,:);
        DEV.ER = DEV_ER(layer,:);
        DAT = tmm1d(DEV,SRC);
        REF(layer, nlam) = DAT.REF;
        TRN(layer, nlam) = DAT.TRN;
        CON(layer, nlam) = DAT.CON;

        REFdB(layer,nlam) = DAT.RdB;
        TRNdB(layer,nlam) = DAT.TdB;
        CONdB(layer,nlam) = DAT.CdB;
    end
    %plotLinear(REF(layer,:), TRN(layer,:), CON(layer,:), LAMBDA, NLAM);
end

%Binary Processing
%REFBinary = 1*(REF >= 0.5) + 0*(REF < 0.5);
%TRNBinary = 1*(TRN >= 0.5) + 0*(TRN < 0.5);

DATAREF = zeros(length(DEVThick(:,1))*length(REF(1,:)),2+size(DEVThick, 2)+size(DEV_ER, 2));
DATATRN = zeros(length(DEVThick(:,1))*length(TRN(1,:)),2+size(DEVThick, 2)+size(DEV_ER, 2));
factor = 0;
for j = 1:length(DEVThick(:,1))
    for i = 1:length(REF(1,:))
        DATAREF(i+factor*100,:) = [LAMBDA(i) DEVThick(j,:) DEV_ER(j,:) REF(j, i)];
        DATATRN(i+factor*100,:) = [LAMBDA(i) DEVThick(j,:) DEV_ER(j,:) TRN(j, i)];
    end
    factor = factor + 1;
end

DATAREF_wavelength1 = DATAREF(1:100:end,:);
DATATRN_wavelength1 = DATATRN(1:100:end,:);

DATAREF_WX_thick = DATAREF_wavelength1(:,2:5)./3;
DATAREF_WX_ER = DATAREF_wavelength1(:,6:9)./25;
DATAREF_WX = round([DATAREF_WX_thick DATAREF_WX_ER],3);
DATAREF_WY = DATAREF_wavelength1(:,10);

DATATRN_WX_thick = DATATRN_wavelength1(:,2:5)./3;
DATATRN_WX_ER = DATATRN_wavelength1(:,6:9)./25;
DATATRN_WX = round([DATATRN_WX_thick DATATRN_WX_ER],3);
DATATRN_WY = DATATRN_wavelength1(:,10);

%Full Dataset, with frequency

DATAREF_X = DATAREF(:,1:9);
DATAREF_Y = DATAREF(:,10);

%Processing for Binary Classification:
%For Wavelength = 1
DATAREF_WYB = round(DATAREF_WY*10); 
logistic_range = linspace(1,10,10);
DATAREF_WYB_Logistic = zeros(size(DATAREF_WYB,1),length(logistic_range));

for index = 1:size(DATAREF_WYB,1)
    DATAREF_WYB_Logistic(index,:) = (DATAREF_WYB(index,:) == logistic_range);
end

%For all wavelengths:
DATAREF_YB = round(DATAREF_Y*10); 
logistic_range_full = linspace(1,10,10);
DATAREF_YB_Logistic = zeros(size(DATAREF_YB,1),length(logistic_range_full));

for index = 1:size(DATAREF_WYB,1)
    DATAREF_YB_Logistic(index,:) = (DATAREF_YB(index,:) == logistic_range_full);
end


save("DATAREF"+num2str(seed)+"RAWFile.csv", "DATAREF", "-ascii");
save("DATAREF"+num2str(seed)+"W1File.csv", "DATAREF_wavelength1", "-ascii");

save("DATATRN"+num2str(seed)+"RAWFile.csv", "DATATRN", "-ascii");
save("DATATRN"+num2str(seed)+"W1File.csv", "DATATRN_wavelength1", "-ascii");

save("DATAREF_WX"+num2str(seed)+"W1File.csv", "DATAREF_WX", "-ascii");
save("DATAREF_WY"+num2str(seed)+"W1File.csv", "DATAREF_WY", "-ascii");


%% Neural Network
% Solve a Pattern Recognition Problem with a Neural Network
% Script generated by Neural Pattern Recognition app
% Created 08-Apr-2021 02:22:33
%
% This script assumes these variables are defined:
%
%   DATAREF_WX - input data.
%   DATAREF_WYB_Logistic - target data.

x = DATAREF_WX';
t = DATAREF_WYB_Logistic';

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.

% Create a Pattern Recognition Network
hiddenLayerSize = 10000;
net = patternnet(hiddenLayerSize, trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 80/100;
net.divideParam.valRatio = 10/100;
net.divideParam.testRatio = 10/100;

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y) %(mean-square-error) 
tind = vec2ind(t);
yind = vec2ind(y);
percentErrors = 100*sum(tind ~= yind)/numel(tind)

% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotconfusion(t,y)
%figure, plotroc(t,y)


[~,ysim] = max(net(DATAREF_WX'));
Accuracy = 100*sum(DATAREF_WYB == ysim')/length(DATAREF_WYB);

%Try a single sample:
sampleDevice = [[1.9 2.9 2.5 1.6]/3 [1.5.^2 2.1.^2 1.9.^2 4.1.^2]/25];
my_sampleDevice_label = 4;
[~,predicted] = max(net(sampleDevice.'));
sampleAccuracy = 100*sum(my_sampleDevice_label == predicted')/length(my_sampleDevice_label);

disp("The training accuracy of the Neural Network in predicting Reflection is: " +num2str(Accuracy)+"%");
disp("The sample device accuracy of the Neural Network in predicting Reflection is: " +num2str(sampleAccuracy)+"%");
disp("The sample's true value: " +num2str(my_sampleDevice_label));
disp("The predicted value: " +num2str(predicted'));

%% PLOTTING DATA

%plotLinear(REF, TRN, CON, LAMBDA, NLAM);

%plotdB(REFdB, TRNdB, CONdB, LAMBDA, NLAM);

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
    
    DAT.RdB = 10*log10(DAT.REF);
    DAT.TdB = 10*log10(DAT.TRN);
    DAT.CdB = 10*log10(DAT.CON);
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

function plotL = plotLinear(ref, trn, con, xAxis, numEL)
    f1 = ref;
    f2 = trn;
    f3 = con;
    lineWidth = 1.5;
    % OPEN FIGURE WINDOW
    F1 = figure('Color','w', 'Position', [1922 42 958 954]);
    set(F1, 'Name', 'Wavelength Response Graph');
    set(F1, 'NumberTitle', 'off');
    % PLOT DATA
    %subplot(2,2,1);
    h = plot(xAxis,f1,'-r','LineWidth',lineWidth);
    hold on;
    %subplot(2,2,2);
    plot(xAxis,f2,'--b','LineWidth',lineWidth);
    hold on;
    %subplot(2,2,[3:4]);
    plot(xAxis,f3,'-.k','LineWidth',lineWidth);
    hold off;
    % SET AXIS LIMITS
    xlim([1 xAxis(numEL)]);
    ylim([0 ,1.02]);
    % MAKE LINES THICK AND FONTS BIGGER
    h2 = get(h,'Parent');
    set(h2,'LineWidth',3,'FontSize',18);
    % SET TICK MARKS

    xm = [0:0.1:xAxis(numEL)];
    xt = {};
    for m = 1 : length(xm)
    xt{m} = num2str(xm(m),'%3.1f');
    end
    set(h2,'XTick',xm,'XTickLabel',xt);
    ym = [0.1:0.1:+1];
    yt = {};
    for m = 1 : length(ym)
    yt{m} = num2str(ym(m),'%2.1f');
    end
    set(h2,'YTick',ym,'YTickLabel',yt);

    % LABEL AXES
    xlabel('$ \textrm{Free Space Wavelength}(\mu m) $','Interpreter','latex');
    ylabel('$ \textrm{Response(Linear)} $','Interpreter','latex');
    title('$ \textrm{Plot of the Wavelength Response of the Bragg Grating}$',...
        'Interpreter', 'LaTeX', 'FontSize', 18);
    % ADD LEGEND
    h = legend('Reflectance','Transmittance', 'Conservation','Location','NorthEastOutside');
    set(h,'LineWidth',2);
    
    % LABEL MINIMUM
    %text(1.5, 0.95, 'Reflectance', 'Color', 'b', 'HorizontalAlignment', 'center');
    %text(1.5, 0.05, 'Transmittance', 'Color', 'r', 'HorizontalAlignment', 'center');
end

function plotdB1 = plotdB(ref, trn, con, xAxis, numEL)
    f1 = ref;
    f2 = trn;
    f3 = con;
    lineWidth = 1.5;
    % OPEN FIGURE WINDOW
    F1 = figure('Color','w', 'Position', [2882 42 958 954]);
    set(F1, 'Name', 'Wavelength Response Graph');
    set(F1, 'NumberTitle', 'off');
    % PLOT DATA
    %subplot(2,2,1);
    h = plot(xAxis,f1,'-r','LineWidth',lineWidth);
    hold on;
    %subplot(2,2,2);
    plot(xAxis,f2,'--b','LineWidth',lineWidth);
    hold on;
    %subplot(2,2,[3:4]);
    plot(xAxis,f3,'-.k','LineWidth',lineWidth);
    hold off;
    % SET AXIS LIMITS
    xlim([1 xAxis(numEL)]);
    ylim([-30 0.5]);
    % MAKE LINES THICK AND FONTS BIGGER
    h2 = get(h,'Parent');
    set(h2,'LineWidth',3,'FontSize',18);
    % SET TICK MARKS

    xm = [0:0.1:xAxis(numEL)];
    xt = {};
    for m = 1 : length(xm)
    xt{m} = num2str(xm(m),'%3.1f');
    end
    set(h2,'XTick',xm,'XTickLabel',xt);
    ym = [-30:5:0];
    yt = {};
    for m = 1 : length(ym)
    yt{m} = num2str(ym(m),'%2.1f');
    end
    set(h2,'YTick',ym,'YTickLabel',yt);

    % LABEL AXES
    xlabel('$ \textrm{Free Space Wavelength}(\mu m) $','Interpreter','latex');
    ylabel('$ \textrm{Response(dB)} $','Interpreter','latex');
    title('$ \textrm{Plot of the Wavelength Response of the Bragg Grating}$',...
        'Interpreter', 'LaTeX', 'FontSize', 18);
    % ADD LEGEND
    h = legend('Reflectance','Transmittance', 'Conservation','Location','NorthEastOutside');
    set(h,'LineWidth',2);
    
    % LABEL MINIMUM
    %text(1.5, 0.95, 'Reflectance', 'Color', 'b', 'HorizontalAlignment', 'center');
    %text(1.5, 0.05, 'Transmittance', 'Color', 'r', 'HorizontalAlignment', 'center');
end

