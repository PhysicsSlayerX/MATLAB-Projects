function DAT = rcwa3d(DEV, SRC)
% RCWA3D            Three-Dimensional Rigorous Coupled-Wave Analysis
%
% DAT = rcwa3d(DEV, SRC);
%
% INPUT ARGUMENTS
% ================
% DEV.er1           relative permittivity in reflection region
% DEV.ur1           relative permeability in reflection region
% DEV.er2           relative permittivity in transmission region
% DEV.ur2           relative permeability in transmission region
%
% DEV.t1            first direct lattice vector
% DEV.t2            second direct lattice vector
%
% DEV.ER            3D array containing convolution matrices for permittivity
%                   N1 by N2 by NLAY, where NLAY is the number of layers.
%
% DEV.UR            3D array containing convolution matrices for permeability
%                   N1 by N2 by NLAY, where NLAY is the number of layers.
%
% DEV.NP            number of spatial harmonics along T1
% DEV.NQ            number of spatial harmonics along T2
% DEV.L             array containing thickness of all layers
%
% SRC.lam()         free space wavelength
% SRC.theta         elevation angle of incidence
% SRC.phi           azimuthal angle of incidence
% SRC.pte           complex amplitude of TE polarization
% SRC.pte           complex amplitude of TM polarization

% OUTPUT ARGUMENTS
% ================
% DAT.RDE           array of diffraction efficiencies of reflected waves
% DAT.TDE           array of diffraction efficiencies of transmitted waves
% DAT.REF           overall reflectance
% DAT.TRN           overall transmittance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM RCWA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE REFRACTIVE INDEX OF EXTERNAL MEDIUMS
n1 = sqrt(DEV.ur1*DEV.er1);
n2 = sqrt(DEV.ur2*DEV.er2);

% CALCULATE RECIPROCAL LATTICE VECTORS
d  = DEV.t1(1)*DEV.t2(2) - DEV.t2(1)*DEV.t1(2);
T1 = 2*pi*[+DEV.t2(2)/d ; -DEV.t2(1)/d ];
T2 = 2*pi*[-DEV.t1(2)/d ; +DEV.t1(1)/d ];

% CALCULATE WAVE VECTOR EXPANSION
k0     = 2*pi/SRC.lam0;
kinc   = n1*[ sin(SRC.theta)*cos(SRC.phi) ...
            ; sin(SRC.theta)*sin(SRC.phi) ...
            ; cos(SRC.theta) ];
        
p      = [-floor(DEV.NP/2):+floor(DEV.NP/2)];
q      = [-floor(DEV.NQ/2):+floor(DEV.NQ/2)];
[Q, P] = meshgrid(q,p);
Kx     = kinc(1) - P*T1(1)/k0 - Q*T2(1)/k0;
Ky     = kinc(2) - P*T1(2)/k0 - Q*T2(2)/k0;
Kzref  = conj(sqrt(DEV.ur1*DEV.er1 - Kx.^2 - Ky.^2));
Kztrn  = conj(sqrt(DEV.ur2*DEV.er2 - Kx.^2 - Ky.^2));

% FORM DIAGONAL K MATRICES
Kx     = diag(Kx(:));
Ky     = diag(Ky(:));
Kzref  = diag(Kzref(:));
Kztrn  = diag(Kztrn(:));

% BUILD SPECIAL MATRICES
NH = DEV.NP*DEV.NQ;
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
NLAY = length(DEV.L);
for nlay = 1 : NLAY
    
    % Build Eigen-Value Problem
    ur = DEV.UR(:,:,nlay);
    er = DEV.ER(:,:,nlay);
    P  = [ Kx/er*Ky , ur-Kx/er*Kx ; Ky/er*Ky-ur , -Ky/er*Kx ];
    Q  = [ Kx/ur*Ky , er-Kx/ur*Kx ; Ky/ur*Ky-er , -Ky/ur*Kx ];
        
    % Compute Eigen-Modes
    [W,LAM]  = eig(P*Q);
    LAM      = sqrt(LAM);
    V        = Q*W/LAM;
    X        = expm(-LAM*k0*DEV.L(nlay));
    
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
    Q    = (1/DEV.ur1) * [Kx*Ky, DEV.ur1*DEV.er1*I - Kx^2 ...
                         ;Ky^2 - DEV.ur1*DEV.er1*I, -Ky*Kx ];
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
    Q    = (1/DEV.ur2) * [Kx*Ky, DEV.ur2*DEV.er2*I - Kx^2 ...
                         ;Ky^2 - DEV.ur2*DEV.er2*I, -Ky*Kx ];
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
if abs(SRC.theta) < 1e-3
    ate = [0;1;0];
else
    ate = cross(kinc, n);
    ate = ate/norm(ate);
end

atm = cross(ate, kinc);
atm = atm/norm(atm);

EP  = SRC.pte*ate + SRC.ptm*atm;
EP  = EP/norm(EP);

% CALCULATE ELECTRIC FIELD SOURCE VECTOR
delta     = zeros(NH, 1);
p0        = ceil(DEV.NP/2);
q0        = ceil(DEV.NQ/2);
m0        = (q0 - 1)*DEV.NP + p0;
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
DAT.RDE       = abs(rx).^2 + abs(ry).^2 + abs(rz).^2;
DAT.RDE       = real(Kzref/kinc(3))*DAT.RDE;
DAT.RDE       = reshape(DAT.RDE, DEV.NP, DEV.NQ);

DAT.TDE       = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;
DAT.TDE       = real(DEV.ur1/DEV.ur2*Kztrn/kinc(3))*DAT.TDE;
DAT.TDE       = reshape(DAT.TDE, DEV.NP, DEV.NQ);

% CALCULATE OVERALL REFLECTANCE AND TRANSMITTANCE
DAT.REF       = sum(DAT.RDE(:));
DAT.TRN       = sum(DAT.TDE(:));


end







