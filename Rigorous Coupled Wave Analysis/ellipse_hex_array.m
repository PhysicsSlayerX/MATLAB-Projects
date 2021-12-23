% runsim_conicaldiffraction.m

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

% SET CONVERGENCE STUDY
convergence_study = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
SRC.lam0  = 1000 * nanometers;
SRC.theta = 30*degrees;
SRC.phi   = 45*degrees;
SRC.pte   = 1/sqrt(2);
SRC.ptm   = 1i/sqrt(2);

% DEVICE PARAMETERS
n_air  = 1.0;
n_grat = 1000;
n_fs   = 1.5100;
a      = 2.5 * micrometers;
r      = 0.50*a;
d      = 500 * nanometers;


DEV.er1 = 1.0;
DEV.ur1 = 1.0;
DEV.er2 = n_grat^2;
DEV.ur2 = 1.0;

DEV.t1 = [ a/2 ; -a*sqrt(3)/2 ];
DEV.t2 = [ a/2 ; +a*sqrt(3)/2 ];

% RCWA PARAMETERS
N1     = 512;
N2     = N1;
DEV.NP = 21;
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

ax     = 0.5; 
ay     = 1;

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
%ER      = 1 - ER;
ERR     = n_air^2 + (n_grat^2 - n_air^2)*ER;
URR     = ones(N1, N2);
DEV.L   = d;

% SHOW DEVICE
subplot(1,1,1);
pcolor(XO, YO, ERR);
shading interp;
axis equal tight;
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONVERGENCE STUDY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if convergence_study
    
    % DEFINE CONVERGENCE STUDY PARAMETERS
    PQ_DAT = [1 : 2 : 31];
    NDAT   = length(PQ_DAT);
    
    % INITIALIZE DATA RECORDS
    REF = zeros(1, NDAT);
    TRN = zeros(1, NDAT);
    
    % 
    % MAIN LOOP
    %
    
    for ndat = 1 : NDAT
        
        % Get Next P & Q
        DEV.NP = PQ_DAT(ndat);
        DEV.NQ = PQ_DAT(ndat);
        
        % COMPUTE CONVOLUTION MATRICES
        NLAY     = length(DEV.L);
        NH       = DEV.NP*DEV.NQ;
        DEV.ER   = zeros(NH,NH, NLAY);
        DEV.UR   = zeros(NH,NH, NLAY);

        for nlay   = 1 : NLAY
            DEV.ER(:,:,nlay) = convmat(ERR(:,:,nlay),DEV.NP,DEV.NQ);
            DEV.UR(:,:,nlay) = convmat(URR(:,:,nlay),DEV.NP,DEV.NQ);
        end
        
        % CALL RCWA3D
        DAT = rcwa3d(DEV, SRC);
        
        % Record Results
        REF(ndat) = DAT.REF;
        TRN(ndat) = DAT.TRN;
        CON       = REF + TRN;
        
        % Show Results
        plot(PQ_DAT(1:ndat), REF(1:ndat), '-or');
        xlim([PQ_DAT(1)-1 , PQ_DAT(NDAT)+1]);
        ylim([0.1 0.3]);
        title('Convergence Study');
        
        drawnow;
    end
    
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALL RCWA3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% COMPUTE CONVOLUTION MATRICES
NLAY = length(DEV.L);
NH   = DEV.NP*DEV.NQ;
DEV.ER   = zeros(NH,NH, NLAY);
DEV.UR   = zeros(NH,NH, NLAY);

for nlay   = 1 : NLAY
    DEV.ER(:,:,nlay) = convmat(ERR(:,:,nlay),DEV.NP,DEV.NQ);
    DEV.UR(:,:,nlay) = convmat(URR(:,:,nlay),DEV.NP,DEV.NQ);
end

% CALL FUNCTION
DAT = rcwa3d(DEV, SRC);

CON = DAT.REF + DAT.TRN;

% REPORT RESULTS
disp(['REF = ' num2str(100*DAT.REF, '%6.2f') '%']);
disp(['TRN = ' num2str(100*DAT.TRN, '%6.2f') '%']);
disp('================');
disp(['CON = ' num2str(100*CON, '%6.2f') '%']);


% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);







