% rcwa3d_proto.m

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
lam0  = 1540 * nanometers;
theta = 0*degrees;
phi   = 0*degrees;
pte   = 1;
ptm   = 0;

% DEVICE PARAMETERS
n_SiO = 1.4496;
n_SiN = 1.9360;
n_fs  = 1.5100;
a     = 1150 * nanometers;
r     =  400 * nanometers;
h1    =  230 * nanometers;
h2    =  345 * nanometers;

er1 = 1.0;
ur1 = 1.0;
er2 = n_fs^2;
ur2 = 1.0;

t1 = [ a/2 ; -a*sqrt(3)/2 ];
t2 = [ a/2 ; +a*sqrt(3)/2 ];

% RCWA PARAMETERS
N1 = 512;
N2 = N1;
NP = 21;
NQ = NP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBLIQUE MESHGRID
p      = linspace(-0.5, +0.5, N1);
q      = linspace(-0.5, +0.5, N2);
[Q, P] = meshgrid(q, p);
XO     = P*t1(1) + Q*t2(1);
YO     = P*t1(2) + Q*t2(2);

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
L   = h1;

% ADD ADDITIONAL LAYERS
ERR(:,:,2) = (n_SiN^2)*ones(N1, N2);
URR(:,:,2) = ones(N1, N2);
L          = [L h2];

% COMPUTE CONVOLUTION MATRICES
NLAY = length(L);
NH   = NP*NQ;
ER   = zeros(NH,NH, NLAY);
UR   = zeros(NH,NH, NLAY);

for nlay   = 1 : NLAY
    ER(:,:,nlay) = convmat(ERR(:,:,nlay),NP,NQ);
    UR(:,:,nlay) = convmat(URR(:,:,nlay),NP,NQ);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM RCWA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE REFRACTIVE INDEX OF EXTERNAL MEDIUMS
n1 = sqrt(ur1*er1);
n2 = sqrt(ur2*er2);

% CALCULATE RECIPROCAL LATTICE VECTORS
d  = t1(1)*t2(2) - t2(1)*t1(2);
T1 = 2*pi*[+t2(2)/d ; -t2(1)/d ];
T2 = 2*pi*[-t1(2)/d ; +t1(1)/d ];

% CALCULATE WAVE VECTOR EXPANSION
k0     = 2*pi/lam0;
kinc   = n1*[ sin(theta)*cos(phi); sin(theta)*sin(phi) ; cos(theta) ];
p      = [-floor(NP/2):+floor(NP/2)];
q      = [-floor(NQ/2):+floor(NQ/2)];
[Q, P] = meshgrid(q,p);
Kx     = kinc(1) - P*T1(1)/k0 - Q*T2(1)/k0;
Ky     = kinc(2) - P*T1(2)/k0 - Q*T2(2)/k0;
Kzref  = conj(sqrt(ur1*er1 - Kx.^2 - Ky.^2));
Kztrn  = conj(sqrt(ur2*er2 - Kx.^2 - Ky.^2));

% FORM DIAGONAL K MATRICES
Kx     = diag(Kx(:));
Ky     = diag(Ky(:));
Kzref  = diag(Kzref(:));
Kztrn  = diag(Kztrn(:));

% BUILD SPECIAL MATRICES
I = eye(NH, NH);
Z = zeros(NH, NH);

% CALCULATE EIGEN-MODES OF THE GAP MEDIUM
Kz  = conj(sqrt(I - Kx^2 - Ky^2));
Q   = [Kx*Ky, I-Kx^2 ; Ky^2-I, -Kx*Ky];
W0  = [ I Z ; Z I ];
LAM = [ 1i*Kz Z ; Z 1i*Kz ];

V0  = Q/LAM;

%  INITIALIZE GLOBAL SCATTERING MATRIX
SG.S11 = zeros(2*NH, 2*NH);
SG.S12 = eye(2*NH, 2*NH);
SG.S21 = eye(2*NH, 2*NH);
SG.S22 = zeros(2*NH, 2*NH);

%
% MAIN LOOP -- ITERATE THROUGH THE LAYERS
%
for nlay = 1 : NLAY
    
    % Build Eigen-Value Problem
    ur = UR(:,:,nlay);
    er = ER(:,:,nlay);
    P  = [ Kx/er*Ky , ur-Kx/er*Kx ; Ky/er*Ky-ur , -Ky/er*Kx ];
    Q  = [ Kx/ur*Ky , er-Kx/ur*Kx ; Ky/ur*Ky-er , -Ky/ur*Kx ];
        
    % Compute Eigen-Modes
    [W,LAM]  = eig(P*Q);
    LAM      = sqrt(LAM);
    V        = Q*W/LAM;
    X        = expm(-LAM*k0*L(nlay));
    
    % Calculate Layer Scattering Matrix
    A     = W\W0 + V\V0;
    B     = W\W0 - V\V0;
    D     = A - X*B/A*X*B;
    S.S11 = D\(X*B/A*X*A - B);
    S.S12 = D\X*(A - B/A*B);
    S.S21 = S.S12;
    S.S22 = S.S11;
    
    % Update Global Scattering Matrix
    SG = star(SG, S);
end

% CONNECT TO REFLECTION REGION
    
    % Calculate the Eigen-Modes
    Q    = (1/ur1) * [Kx*Ky, ur1*er1*I - Kx^2 ...
                     ;Ky^2 - ur1*er1*I, -Ky*Kx ];
    Wref = [I Z ; Z I];
    LAM  = [1i*Kzref Z ; Z 1i*Kzref ];
    Vref = Q/LAM;
    
    % Calculate Reflection-Side Scattering Matrix
    A     = W0\Wref + V0\Vref;
    B     = W0\Wref - V0\Vref;
    
    S.S11 = -A\B;
    S.S12 = 2*inv(A);
    S.S21 = 0.5*(A - B/A*B);
    S.S22 = B/A;
    
    % Update Global Scattering Matrix
    SG = star(S, SG);
    
% CONNECT TO TRANSMISSION REGION
    
    % Calculate the Eigen-Modes
    Q    = (1/ur2) * [Kx*Ky, ur2*er2*I - Kx^2 ...
                     ;Ky^2 - ur2*er2*I, -Ky*Kx ];
    Wtrn = [I Z ; Z I];
    LAM  = [1i*Kztrn Z ; Z 1i*Kztrn ];
    Vtrn = Q/LAM;
    
    % Calculate Transmission-Side Scattering Matrix
    A     = W0\Wtrn + V0\Vtrn;
    B     = W0\Wtrn - V0\Vtrn;
    
    S.S11 = B/A;
    S.S12 = 0.5*(A - B/A*B);
    S.S21 = 2*inv(A);
    S.S22 = -A\B;
    
    % Update Global Scattering Matrix
    SG = star(SG, S);
    
% COMPUTE POLARIZATION VECTOR
n = [0;0;1];
if abs(theta) < 1e-3
    ate = [0;1;0];
else
    ate = cross(kinc, n);
    ate = ate/norm(ate);
end

atm = cross(ate, kinc);
atm = atm/norm(atm);

EP  = pte*ate + ptm*atm;
EP  = EP/norm(EP);

% CALCULATE ELECTRIC FIELD SOURCE VECTOR
delta     = zeros(NH, 1);
p0        = ceil(NP/2);
q0        = ceil(NQ/2);
m0        = (q0 - 1)*NP + p0;
delta(m0) = 1;
esrc      = [ EP(1)*delta ; EP(2)*delta ];

% CALCULATE SOURCE VECTORS
csrc      = Wref\esrc;

% CALCULATE REFLECTED FIELDS
cref      = SG.S11*csrc;
eref      = Wref*cref;
rx        = eref(1:NH);
ry        = eref(NH+1 : 2*NH);
rz        = -Kx/Kzref*rx - Ky/Kzref*ry;


% CALCULATE TRANSMITTED FIELDS
ctrn      = SG.S21*csrc;
etrn      = Wtrn*ctrn;
tx        = etrn(1:NH);
ty        = etrn(NH+1 : 2*NH);
tz        = -Kx/Kztrn*tx - Ky/Kztrn*ty;

% CALCULATE DIFFRACTION EFFICIENCIES
RDE       = abs(rx).^2 + abs(ry).^2 + abs(rz).^2;
RDE       = real(Kzref/kinc(3))*RDE;
RDE       = reshape(RDE, NP, NQ);

TDE       = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;
TDE       = real(ur1/ur2*Kztrn/kinc(3))*TDE;
TDE       = reshape(TDE, NP, NQ);

% CALCULATE OVERALL REFLECTANCE AND TRANSMITTANCE
REF       = sum(RDE(:));
TRN       = sum(TDE(:));
CON       = REF + TRN;


% REPORT RESULTS
disp(['REF = ' num2str(100*REF, '%6.2f') '%']);
disp(['TRN = ' num2str(100*TRN, '%6.2f') '%']);
disp('================');
disp(['CON = ' num2str(100*CON, '%6.2f') '%']);









