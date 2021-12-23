% pwem2d_bands.m
% EMPossible

% INITIALIZE MATLAB
close all;
clc;
clear all; 

time1 = clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATTICE PARAMETERS
a      = 1;
erfill = 12.0;
erhole = 1.0;
rhole  = 0.48*a;

% PWEM PARAMETERS
N1 = 1024;
N2 = N1;
NP = 11;
NQ = NP;

% BAND DIAGRAM PARAMETERS
wn_max = 1;
NBETA  = 50;

% ACTIVE PLOT CONTROL
active_plot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTRUCT BETA PATH AROUND IBZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATTICE VECTORS
t1 = a * [0.5 ; -0.5*sqrt(3)];
t2 = a * [0.5 ; +0.5*sqrt(3)];
T1 = (2*pi/a) * [1 ; -1/sqrt(3)];
T2 = (2*pi/a) * [1 ; +1/sqrt(3)];

% DEFINE KEY POINTS OF SYMMETRY
G = [0 ; 0];
M = 0.5*T2;
K = (1/3)*T1 + (1/3)*T2;

% DEFINE ORDER OF KEY POINTS
KP = [G M K G];                   % Key points  
KL = {'\Gamma' 'M' 'K' '\Gamma'}; % Key labels

% CALCULATE BETA AXIS RESOLUTION
L = 0;
NKP = length(KP(1, :));

for nkp = 1 : NKP-1
    L = L + norm(KP(:, nkp + 1) - KP(:, nkp));
end
res = L/NBETA;

% BUILD BETA LIST
BETA = KP(:,1);
KT   = 1;

for nkp = 1 : NKP - 1
    % Get End Points
    kp1 = KP(:, nkp);
    kp2 = KP(:, nkp+1);
    
    % Calculate Number of Points
    L = norm(kp2 - kp1);
    NB = round(L/res);
    
    % Calculate Points from kp1 to kp2
    bx = kp1(1) + (kp2(1) - kp1(1))*[1:NB]/NB;
    by = kp1(2) + (kp2(2) - kp1(2))*[1:NB]/NB;
    
    % Append Points to BETA
    BETA = [ BETA, [bx; by] ];
    KT(nkp + 1) = length(BETA(1,:));
    
end
NBETA = length(BETA(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL AND CONVOLUTION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBLIQUE MESHGRID
p = linspace(-0.5, +0.5, N1);
q = linspace(-0.5, +0.5, N2);
[Q, P] = meshgrid(q, p);
XO = P*t1(1) + Q*t2(1);
YO = P*t1(2) + Q*t2(2);

% BUILD UNIT CELL
b  = a*sqrt(3);
ER = zeros(N1, N2);
ER = ER | ((XO-a/2).^2 + YO.^2) <= rhole^2; 
ER = ER | ((XO+a/2).^2 + YO.^2) <= (rhole/2)^2; 
ER = ER | (XO.^2 + (YO-b/2).^2) <= (rhole/3)^2; 
ER = ER | (XO.^2 + (YO+b/2).^2) <= (rhole/4)^2; 

% CONVERT TO REAL MATERIALS
ER = erfill + (erhole - erfill)*ER;

% CONSTRUCT CONVOLUTIOM MATRIX
ERC = convmat(ER, NP, NQ);

% SHOW DATA
clf;
subplot(131);
pcolor(XO, YO, double(ER));
shading interp;
axis equal tight;
colorbar;
title('UNIT CELL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PHOTONIC BANDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE DATA RECORD
M     = NP*NQ;
WN_TE = zeros(M, NBETA);
WN_TM = zeros(M, NBETA);

% FOURIER COEFFICIENT INDEX MESHGRID
p = [-floor(NP/2): +floor(NP/2)];
q = [-floor(NQ/2): +floor(NQ/2)];
[Q, P] = meshgrid(q,p);

% 
% MAIN LOOP -- ITERATE OVER LIST OF BETA's
%

for nb = 1 : NBETA
    
    % Build K Matrices
    Kx = BETA(1, nb) - P*T1(1) - Q*T2(1);
    Ky = BETA(2, nb) - P*T1(2) - Q*T2(2);
    Kx = diag(Kx(:));
    Ky = diag(Ky(:));
    
    % Calculate TM Bands
    A  = Kx^2 + Ky^2;
    k0 = eig(A, ERC);
    k0 = sort(real(sqrt(k0)));
    WN_TM(:, nb) = a*k0/(2*pi);
    
    % Calculate TE Bands
    A  = Kx/ERC*Kx + Ky/ERC*Ky;
    k0 = eig(A);
    k0 = sort(real(sqrt(k0)));
    WN_TE(:, nb) = a*k0/(2*pi);
    
    if active_plot
        % SUBPLOT
        subplot(1,3, [2 3]);

        % DRAW BANDS
        plot([1:NBETA], WN_TE, '-b','Linewidth', 1); hold on;
        plot([1:NBETA], WN_TM, '-r','Linewidth', 1); hold off;

        % BETA AXIS
        set(gca, 'XTick', KT, 'XTickLabel', KL);
        xlabel('Block Wave Vector $\vec{\beta}$', 'Interpreter', 'LaTex', 'FontSize', 14);

        % WN_AXIS
        ylabel('$\frac{a}{\lambda}$', 'Interpreter', 'LaTex',...
           'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 16);

        % SET GRAPHICS VIEW
        axis tight;
        ylim([0 wn_max]);
        title('PHOTONIC BAND DIAGRAM', 'Interpreter', 'LaTeX');

        drawnow;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE A PROFESSIONAL LOOKING BAND DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SUBPLOT
subplot(1,3, [2 3]);

% DRAW BANDS
plot([1:NBETA], WN_TE, '-b','Linewidth', 1); hold on;
plot([1:NBETA], WN_TM, '-r','Linewidth', 1); hold off;

% BETA AXIS
set(gca, 'XTick', KT, 'XTickLabel', KL);
xlabel('Block Wave Vector $\vec{\beta}$', 'Interpreter', 'LaTex', 'FontSize', 14);

% WN_AXIS
ylabel('$\frac{a}{\lambda}$', 'Interpreter', 'LaTex',...
       'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 16);

% SET GRAPHICS VIEW
axis tight;
ylim([0 wn_max]);
title('PHOTONIC BAND DIAGRAM', 'Interpreter', 'LaTeX');

time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);



