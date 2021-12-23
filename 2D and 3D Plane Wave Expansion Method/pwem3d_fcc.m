% pwem3d_fcc.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees = pi/180;

% OPEN FIGURE WINDOW
fig = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0.05 1 0.95]);
set(fig, 'Name', 'PWEM-3D');
set(fig, 'NumberTitle', 'off');

% SHOW SINGLE OR LATTICE
single_cell = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FCC PARAMETERS
a  = 1;
rs = 0.38*a;
nL = 1.0;
nH = 3.0;

% PWEM PARAMETERS
N1 = 128;
N2 = N1;
N3 = N2;

% NUMBER OF SPATIAL HARMONICS
NP = 3;
NQ = NP;
NR = NP;

NSOL = 50;

% BAND DIAGRAM PARAMETERS
NBETA  = 50;
wn_max = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PATH AROUND THE IBZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRECT LATTICE VECTORS
t1 = a*[0;0.5;0.5];
t2 = a*[0.5;0;0.5];
t3 = a*[0.5;0.5;0];

% RECIPROCAL LATTICE VECTORS
T1 = (2*pi/a)*[-1 ;  1 ;  1];
T2 = (2*pi/a)*[ 1 ; -1 ;  1];
T3 = (2*pi/a)*[ 1 ;  1 ; -1];

% CALCULATE KEY POINTS OF SYMMETRY
PG = [0 ; 0 ; 0];
PX = 0.5*T2 + 0.5*T3;
PL = 0.5*T1 + 0.5*T2 + 0.5*T3;
PK = (3/8)*T1 + (3/8)*T2 + (3/4)*T3;
PU = (1/4)*T1 + (5/8)*T2 + (5/8)*T3;
PW = (1/4)*T1 + (1/2)*T2 + (3/4)*T3;

% DEFINE PATH AROUND IBZ
KP = [ PU PX PG PL PK PW PU ];
KL = {'U' 'X' '\Gamma' 'L' 'K' 'W' 'U'};

% CALCULATE BETA AXIS RESOLUTION
L = 0;
NKP = length(KP(1,:));
for nkp = 1 : NKP-1
    L = L + norm(KP (:, nkp+1) - KP(:, nkp));
end
res = L/NBETA;

% BUILD BETA LIST
BETA = KP(:,1);
KT   = 1;

for nkp = 1 : NKP-1
    % Get End Points
    kp1 = KP(:, nkp);
    kp2 = KP(:, nkp+1);
    % Calculate Number of Points
    L = norm(kp2 - kp1);
    NB = round(L/res);
    % Calculate Beta Points From kp1 to kp2
    bx = kp1(1) + (kp2(1) - kp1(1))*[1:NB]/NB;
    by = kp1(2) + (kp2(2) - kp1(2))*[1:NB]/NB;
    bz = kp1(3) + (kp2(3) - kp1(3))*[1:NB]/NB;
    % Append Points to BETA
    BETA = [ BETA,  [bx;by;bz] ];
    KT(nkp+1) = length(BETA(1,:));
end
NBETA = length(BETA(1,:));

% CHECK PATH AROUND IBZ
% line(BETA(1, :), BETA(2,:), BETA(3,:), 'LineWidth', 1);
% axis equal tight;
% view(110,20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL AND CONVOLUTION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BUILD OBLIQUE MESHGRID
p = linspace(-0.5, + 0.5, N1);
q = linspace(-0.5, + 0.5, N2);
r = linspace(-0.5, + 0.5, N3);

[Q, P, R] = meshgrid(q,p,r);
XO = P*t1(1) + Q*t2(1) + R*t3(1);
YO = P*t1(2) + Q*t2(2) + R*t3(2);
ZO = P*t1(3) + Q*t2(3) + R*t3(3);

% BUILD UNIT CELL
ER = zeros(N1, N2, N3);

d  = -0.5*t1 - 0.5*t2 - 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;
d  = -0.5*t1 - 0.5*t2 + 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;

d  = -0.5*t1 + 0.5*t2 - 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;
d  = -0.5*t1 + 0.5*t2 + 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;

d  = +0.5*t1 - 0.5*t2 - 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;
d  = +0.5*t1 - 0.5*t2 + 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;

d  = +0.5*t1 + 0.5*t2 - 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;
d  = +0.5*t1 + 0.5*t2 + 0.5*t3;
ER = ER | (XO + d(1)).^2 + (YO + d(2)).^2 + (ZO + d(3)).^2 <= rs^2;

ER = 1 - ER;

% SHOW UNIT CELL

clf;
subplot(1,3,1);
ERS = smooth3(ER);
[F1, V1] = isosurface(XO, YO, ZO, ERS, 0.5);
[F2, V2] = isocaps(XO, YO, ZO, ERS, 0.5);
if single_cell
    h = patch('Faces' , F1, 'Vertices',  V1);
    set(h, 'FaceColor', [0.5 0.7 0.5], 'LineStyle', 'none');
    h = patch('Faces' , F2, 'Vertices',  V2);
    set(h, 'FaceColor', [0.5 0.7 0.5], 'LineStyle', 'none');
else
    for n3 = -0.5 : 0.5
        for n2 = -0.5 : 0.5
            for n1 = -0.5 : 0.5
                v = V1;
                v(:,1) = v(:,1) + n1*t1(1) + n2*t2(1) + n3*t3(1);
                v(:,2) = v(:,2) + n1*t1(2) + n2*t2(2) + n3*t3(2);
                v(:,3) = v(:,3) + n1*t1(3) + n2*t2(3) + n3*t3(3);
                h = patch('Faces' , F1, 'Vertices',  v);
                set(h, 'FaceColor', [0.5 0.7 0.5], 'LineStyle', 'none');
                
                v = V2;
                v(:,1) = v(:,1) + n1*t1(1) + n2*t2(1) + n3*t3(1);
                v(:,2) = v(:,2) + n1*t1(2) + n2*t2(2) + n3*t3(2);
                v(:,3) = v(:,3) + n1*t1(3) + n2*t2(3) + n3*t3(3);
                h = patch('Faces' , F2, 'Vertices',  v);
                set(h, 'FaceColor', [0.5 0.7 0.5], 'LineStyle', 'none');
            end
        end
    end
end

axis equal tight;
view(50, 15);
camlight;
lighting phong;
title('UNIT CELL');

% CONVERT ER TO REAL MATERIALS
ER = nL^2 + (nH^2 - nL^2)*ER;

% BUILD THE CONVOLUTION MATRIX
ERC = convmat(ER, NP, NQ, NR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE PHOTONIC BANDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE DATA RECORDS
NH = NP*NQ*NR;
WN = zeros(NSOL, NBETA);

% FOURIER COEFFICIENTS MESHGRID
p = [-floor(NP/2): +floor(NP/2)];
q = [-floor(NQ/2): +floor(NQ/2)];
r = [-floor(NR/2): +floor(NR/2)];

[Q, P, R] = meshgrid(q, p, r);

%
% MAIN LOOP -- ITERATE OVER BLOCH WAVE VECTOR
%

for nb = 1 : NBETA
    % Calcualte Wave Vector Expansion
    Kx = BETA(1, nb) - P*T1(1) - Q*T2(1) - R*T3(1);
    Ky = BETA(2, nb) - P*T1(2) - Q*T2(2) - R*T3(2);
    Kz = BETA(3, nb) - P*T1(3) - Q*T2(3) - R*T3(3);
    K  = sqrt(abs(Kx).^2 + abs(Ky).^2 + abs(Kz).^2);
    
    % Calculate Perpendicular Polarization Vectors
    p1x = speye(NH, NH);      p2x = speye(NH, NH);
    p1y = speye(NH, NH);      p2y = speye(NH, NH);
    p1z = speye(NH, NH);      p2z = speye(NH, NH);
    for m = 1 : NH
        k = [Kx(m) ; Ky(m) ; Kz(m)];
        if norm(k) < 1e-10
            p1 = [1;0;0];
            p2 = [0;1;0];
        else
            v = [4*k(2) ; 2*k(3) ; 3*k(1)];
            p1 = cross(k,v);
            p1 = p1/norm(p1);
            p2 = cross(k, p1);
            p2 = p2/norm(p2);
        end
        p1x(m,m) = p1(1);
        p1y(m,m) = p1(2);
        p1z(m,m) = p1(3);
        
        p2x(m,m) = p2(1);
        p2y(m,m) = p2(2);
        p2z(m,m) = p2(3);
    end
    
    % Build Eigen-Value Problem
    K = diag(sparse(K(:)));
    A11 = +K*p2x/ERC*K*p2x + K*p2y/ERC*K*p2y + K*p2z/ERC*K*p2z;
    A12 = -K*p2x/ERC*K*p1x - K*p2y/ERC*K*p1y - K*p2z/ERC*K*p1z;
    A21 = -K*p1x/ERC*K*p2x - K*p1y/ERC*K*p2y - K*p1z/ERC*K*p2z;
    A22 = +K*p1x/ERC*K*p1x + K*p1y/ERC*K*p1y + K*p1z/ERC*K*p1z;
    
    A   = [ A11 A12 ; A21 A22 ];
    
    % Solve Eigen-Value Problem
    k02 = sort(eig(A));
    WN(:,nb) = (a/2/pi)*real(sqrt(k02(1: NSOL)));
    
    % Draw Bands
    subplot(1,3, [2 3]);
    plot([1:nb], WN(:, 1 : nb)', 'Color', 'b');
    axis tight;
    xlim([1 NBETA]);
    ylim([0 wn_max]);
    
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW A PROFESSIONAL BAND DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREPARE AXES

cla;
hold on;

% DRAW VERTICAL LINES
for m = 1 : length(KT)
    x = KT(m) * [1 1];
    y = [0 wn_max];
    line(x,y, 'LineStyle', '--', 'Color', 'k');
end

% DRAW BANDS
plot([1:NBETA], WN , 'b', 'LineWidth', 1);

% BETA AXIS
set(gca, 'XTick', KT, 'XTickLabel', KL);
xlabel('Block Wave Vector $\vec{\beta}$', 'Interpreter', 'LaTeX');

% WN AXIS
ylabel('Normalized Frequency $\omega_n = a/\lambda_0$', 'Interpreter', 'LaTeX');
YT = [0:0.2:wn_max];
YL = {'0'};
for m = 2 : length(YT)
    YL{m} = num2str(YT(m), '%3.1f' );
end
set(gca, 'YTick', YT, 'YTickLabel', YL);


% SET GRAPHICS VIEW
hold off;
axis tight;
ylim([0 wn_max]);
title('PHOTONIC BAND DIAGRAM');






