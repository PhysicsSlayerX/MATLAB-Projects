% BASED from the University of Texas - El Paso EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This code was written by Francis S. Dela Cruz, and was heavily referenced
% from the lectures on CEM (Computational Electromagnetics) by Dr. Raymond Rumpf. 
% The code was built from scratch with just the help of CEM Lectures.
% The basis of this code is from the Maxwell's Equation that was translated
% into matrices to be solved in MATLAB. This is a benchmarking code for
% future FDFD Calculations. Note that this was written with just a basic
% knowledge in MATLAB. Optimization and refactoring of codes are necessary
% to keep the runtime lower and simulation speed faster, as the coder
% implements better code syntax and alogrithms.

% Benchmarked as of June 22, 2020

% This MATLAB Program implements the Finite Difference Frequency Domain Method (FDFD).
close all;
clc;
clear all;

% OPEN FIGURE WINDOW
fig = figure('Color','w');
set(gcf, 'Position', [1921 41 960 963]);

% FORMAT SIGNIFICANT DIGITS
format short g;
% UNITS
centimeters = 1;
millimeters = 0.1 * centimeters;
meters = 100 * centimeters;
degrees = pi/180;
seconds = 1;
hertz = 1/seconds;
gigahertz = 1e9 * hertz;
megahertz = 1e6 * hertz;
% CONSTANTS
c0 = 299792458 * meters/seconds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a simulation of a dielectric grating

% SOURCE PARAMETERS
f0 = 24* gigahertz; %operating frequency
lam0 = c0/f0; %operating wavelength in free space
theta = 15* degrees; %angle of incidence
MODE = 'H'; %electromagnetic mode, 'E' or 'H'
% GRATING PARAMETERS
fd = 8.0 * gigahertz; %design frequency
lamd = c0/fd; %design wavelength
x1 = 0.1040*lamd; %width of tooth 1
x2 = 0.0175*lamd; %width of slot
x3 = 0.1080*lamd; %width of tooth 2
L = 0.6755*lamd; %period of grating
d = 0.2405*lamd; %grating depth
t = 0.0510*lamd; %substrate thickness
mech_feat = [x1 x2 x3 L d t]; %array of mechanical features
ur = 1.0; %relative permeability of grating
er = 10.0; %dielectric constant of grating
% EXTERNAL MATERIALS
ur1 = 1.0; %permeability in the reflection region
er1 = 1.0; %permittivity in the reflection region
ur2 = 1.0; %permeability in the transmission region
er2 = 1.0; %permittivity in the transmission region
% GRID PARAMETERS
NRES = 40; %grid resolution
BUFZ = 2*lam0 * [1 1]; %spacer region above and below grating
NPML = [0 0 20 20]; %size of PML at top and bottom of grid
xbc = -2; % xbc tells if the x-axis will have DIrichlet or Pseudo-Periodic/Floquet Boundary Condition
ybc = 0; % xbc tells if the y-axis will have DIrichlet or Pseudo-Periodic/Floquet Boundary Condition

% OPTIMIZED PARAMETERS

% Refractive Indices of the Simulation Setup

n_device = sqrt(ur*er); % refractive index of the device
nref = sqrt(ur1*er1); % refractive index of the reflection region
ntrn = sqrt(ur2*er2); % refractive index of the transmission region

% Consider the Wavelength for resolving the resolution of the simulation

n_max = max([n_device nref ntrn]);
lambda_min = lam0/n_max;
delta_gamma = lambda_min/NRES;

% Consider the Mechanical Features for optimizing the simulaiton

d_min = min(mech_feat);
delta_d = d_min/6;

% Choosing which is smaller between the parameters (to resolve smallest
% features

delta_x = min(delta_gamma, delta_d);
delta_y = delta_x;

% Resolving Critical Dimensions through "Snapping"
xTdev = x1+x2+x3; %the critical dimension here is the duty cycle of the grating
yTdev = d+t; %the height and the substrate thickness

Mx = ceil(xTdev/delta_x);
My = ceil(yTdev/delta_y);

% Total Grid Size

dx = xTdev/Mx; %the adjusted dx
dy = yTdev/My; %the adjusted dy

Nx = round(L/dx);
Nx = 2*round(Nx/2)+1; %physical size along x-axis
Ny = round(yTdev/dy)+ 2.*NPML(3)+ 2.*ceil(2*lam0/(nref*dy)); %physical size along y-axis


%The 2x Grid Parameters

Nx2 = 2*Nx; %The 2x Grid requires 2x of the number of elements for both x and y
Ny2 = 2*Ny;

dx2 = dx/2;
dy2 = dy/2;

% Building the Simulation Matrix and the device

nx1 = floor((Nx2 - round(xTdev/dx2))/2)-1; %location of the left side of the first tooth from the left side of the grid
nx2 = nx1 + round(x1/dx2) - 1; %location of the right side of the first tooth from the left side of the grid
nx3 = nx2+1; %location of the left side of the gap from the left side of the grid
nx4 = nx3 + round(x2/dx2) - 1; %location of the right side of the gap from the left side of the grid
nx5 = nx4+1; %location of the left side of the second tooth from the left side of the grid
nx6 = nx5 + round(x3/dx2) - 1; %location of the right side of the second tooth from the left side of the grid

ny1 = floor((Ny2 - round(yTdev/dy2))/2); %location of the top surface of the grating from the top of the grid
ny2 = ny1 + round(d/dy2) - 1; %location of the tooth base from the top of the grid
ny3 = ny2+1; %location of the top surface of the grating substrate from the top of the grid
ny4 = ny3 + round(t/dy2) - 1; %location of the bottom surface of the grating from the top of the grid

ER2 = er1*ones(Nx2,Ny2); %Initiate a grid of all ones times the dielectric of the reflection region (er1).
ER2(nx1:nx2, ny1:ny2) = er; %Generate and fill the first tooth with the appropriate dielectric values
ER2(nx5:nx6, ny1:ny2) = er; %Generate and fill the second tooth with the appropriate dielectric values
ER2(:,ny3:ny4) = er; %Generate and fill the grating substrate
ER2(:,ny4+1:Ny2) = er2; %Fill the transmission region with its dielectric values

UR2 = ur1*ones(Nx2,Ny2); %Since the relative permeability remains 1 (non-magnetic), fill the entire grid with ur1 = 1

% The Finite Difference Frequency Domain Method

%Wave vector with Source Frequency

k0 = 2*pi./lam0; 

% Material Properties in the Reflected and Transmitted Regions

er_ref = ER2(:, 2*NPML(3) + 10); %One layer in the reflection region's relative permittivity 
er_ref = mean(er_ref(:)); %mean relative permittivity of the grid layer
er_trn = ER2(:, Ny2 - 2*NPML(3) - 5);%One layer in the transmission region's relative permittivity
er_trn = mean(er_trn(:)); %mean relative permittivity of the grid layer

ur_ref = UR2(:, 2*NPML(4) + 10); %One layer in the reflection region's relative permeability 
ur_ref = mean(ur_ref(:)); %mean relative permittivity of the grid layer
ur_trn = UR2(:, Ny2 - 2*NPML(4) - 5);%One layer in the transmission region's relative permeability
ur_trn = mean(ur_trn(:)); %mean relative permittivity of the grid layer

nref(er_ref*ur_ref); %refractive index of the reflection region
ntrn(er_trn*ur_trn); %refractive index of the transmission region

if er_ref <0 && ur_ref<0
    nref = - nref;
end

if er_trn <0 && ur_trn <0
    ntrn = - ntrn;
end
    
% Calculate the Perfectly Matched Layer 

NGRID2X = [Nx2 Ny2];
[sx, sy] = calcpml2d(NGRID2X,2*NPML); % The function that calculates the Perfectly Matched Layer 

% Incorporate the PML to the grid: (Double Diagonally Anisotropic)

URxx = UR2./sx.*sy;
URyy = UR2.*sx./sy;
URzz = UR2.*sx.*sy;
ERxx = ER2./sx.*sy;
ERyy = ER2.*sx./sy;
ERzz = ER2.*sx.*sy;

% Overlay Materials Onto 1x Grids (Based on the Yee Grid)

URxx = URxx(1:2:Nx2,2:2:Ny2);
URyy = URyy(2:2:Nx2,1:2:Ny2);
URzz = URzz(2:2:Nx2,2:2:Ny2);
ERxx = ERxx(2:2:Nx2,1:2:Ny2);
ERyy = ERyy(1:2:Nx2,2:2:Ny2);
ERzz = ERzz(1:2:Nx2,1:2:Ny2);

% Compute the Wave Vector Terms

kinc = k0.*nref.*[sin(theta); cos(theta)]; % Incident Wave Vector

m = [-floor(Nx/2):floor(Nx/2)]';    % Diffraction Modes (Set of Integer)
kx_m = kinc(1) - m*(2*pi/(Nx*dx));  % Transverse Component of the Incident Wave Vector

ky_ref_m = (sqrt((k0*nref)^2 - (kx_m).^2)); % Longitudinal Component, Reflection Region
ky_trn_m = (sqrt((k0*ntrn)^2 - (kx_m).^2)); % Longitudinal Component, Transmission Region

% Convert Diagonal Materials Matrices

% NOTE: USE SPARSE MATRIX, NOT FULL MATRIX!

ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% Construct the Derivative Matrices

NGRID = [Nx Ny]; % 1x Grid Dimensions
RES = [dx dy];   % Grid unit cell's dimensions
BC = [xbc ybc];  % Boundary Condition: 0 for Dirichlet, -2 for Pseudo-Periodic/Floquet

%This function construct derivative matrices for Maxwell's Equation - 
%The Yeeder
[DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0); 

% Compute the Wave Matrix A - the Wave Equation in Matrix Form
% E-Mode and H-Mode

switch MODE
    case 'E'
        A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
    case 'H'
        A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
    otherwise
        error('Unrecognized Polarization');
end

% Compute the Source Field
x1a = m*dx;
y1a = [1:Ny]*dy;
[Y1a, X1a] = meshgrid(y1a, x1a);
f_src = exp(1i.*(kinc(1).*X1a + kinc(2).*Y1a)); %Source Field Matrix
fsrc = f_src(:); %reshaping the Source Field Matrix into a column vector

% Compute the Scattered-Field Masking Matrix, Q

Q = sparse(Nx, Ny); % Initiate the Q Matrix - Scattered-Field Masking Matrix
sr = 1; % Set the amplitude to 1
Q(:, 1:NPML(3)+20) = sr; %Fill the upper part of the 1x Grid with sr values
% leaving the rest of the grid with 0 values. The area with 1's is the
% Masked Area where the scattering fields are. The area with 0's is the 
% total field area. Both areas can be at any arbitrary sizes depending on
% the requirements of your simulation.

Q = diag(sparse(Q(:))); %Initialize the diagonalizing the sparsed Q.

% Compute the source vector, b; 

b = (Q*A - A*Q)*fsrc;

% Calculating the E and H Fields using Af = b, f = A\b (pre-division
% matrix), A\b means inv(A)/b, but don't use inv(A) in this calculation.
% It is not very efficient.

% E-mode - Hz Field
% H-mode - Ez Field
Field = A\b;
Field = full(Field); %shows the full vector, not the sparse one
Field = reshape(Field,Nx,Ny); %reshape the column vector to matrix form.

% POST PROCESSING: Computation of the Diffraction Efficiencies

% Extract Transmitted and Reflected Fields

Fref = Field(:, NPML(3)+10); 
Ftrn = Field(:, Ny - NPML(3)-10);

% Remove the Phase Tilt

% The phase was introduced by the grating. Since the wave will follow the
% phase introduced by the grating's periodicity, we must remove it in the
% reflection and transmission fields. It's like removing the 'tilt' in the
% wave. 

x_value = [1:Nx]'*dx;
phaseTilt = exp(-1i*(kinc(1)*x_value));
Aref = Fref.* phaseTilt;
Atrn = Ftrn.* phaseTilt;

% Calculate the Complex Amplitudes of the Spatial Harmonics

% Spatial Harmonics can be calculated by using the MATLAB's 
% Fast Fourier Transform function in the
% T-R Fields, then use the fftshift to center the values, then flip the
% entire x-row upside down to show the modes. These are the S and U values,
% the Scattering Parameters S11 and S21 for the E Field and U11 and U21 for
% the H Field

Sref = flipud(fftshift(fft(Aref)))/Nx;
Strn = flipud(fftshift(fft(Atrn)))/Nx;

% Calculate Diffraction Efficiencies

% Reflectance and Transmittance per Diffraction Mode

switch MODE
    case 'E'
        R = abs(Sref).^2.*(real(ky_ref_m./ur1)./real(kinc(2)/ur1));
        T = abs(Strn).^2.*(real(ky_trn_m./ur2)./real(kinc(2)/ur1));
    case 'H'
        R = abs(Sref).^2.*(real(ky_ref_m./er1)./real(kinc(2)/er1));
        T = abs(Strn).^2.*(real(ky_trn_m./er2)./real(kinc(2)/er1));
end

% Reflectance and Transmittance - Total of all the respective Diffraction
% Efficiencies

REF = 100.*sum(R);
TRN = 100.*sum(T);

CON = REF + TRN;

% Displaying the Results

index_m = Nx/2 + 1/2;
disp('Reflection Diffraction Orders:' );
disp(['RDE(-1) = '  num2str(100*R(index_m-1))]);
disp(['RDE(0) = '  num2str(100*R(index_m))]);
disp(['RDE(1) = '  num2str(100*R(index_m+1))]);
disp(['RDE(2) = '  num2str(100*R(index_m+2))]);

disp('Transmission Diffraction Orders:' );
disp(['TDE(-1) = '  num2str(100*T(index_m-1))]);
disp(['TDE(0) = '  num2str(100*T(index_m))]);
disp(['TDE(1) = '  num2str(100*T(index_m+1))]);
disp(['TDE(2) = '  num2str(100*T(index_m+2))]);


disp(['REF = ' num2str(REF)]);
disp(['TRN = ' num2str(TRN)]);
disp(['CON = ' num2str(CON)]);


% Plotting the 2D Device
numgraph = 4;

subplot(1,numgraph,1);

xa = [-Nx2/2:Nx2/2]*dx2;
ya = [0:Ny2-1]*dy2;
[Y, X] = meshgrid(ya, xa);

h = imagesc(xa,ya, UR2.',[1 10]);
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\epsilon_r');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;

subplot(1,numgraph,2);

g = imagesc(xa,ya, ER2.');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\mu_r');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;


subplot(1,numgraph,3);
xb = [m]*dx;
yb = [1:Ny]*dy;
[Yb, Xb] = meshgrid(yb, xb);
h = imagesc(xb,yb, real(Field).');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
mode = convertCharsToStrings(MODE);
title([mode+ '-mode at ' + num2str(f0/gigahertz) + ' GHz', 'Re\{F\}']);
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;
%colormap(jet(1024));

subplot(1,numgraph,4);

h = imagesc(xb,yb, imag(Field).');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title([mode+ '-mode at ' + num2str(f0/gigahertz) + ' GHz' , 'Im\{F\}']);
shading interp;
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;
%colormap(jet(1024));


function [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc)
% YEEDER Construct Yee Grid Derivative Operators on a 2D Grid
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc);
%
% Note for normalized grid, use this function as follows:
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);
%
% Input Arguments
% =================
% NGRID [Nx Ny] grid size
% RES [dx dy] grid resolution of the 1X grid
% BC [xbc ybc] boundary conditions
% -2: periodic (requires kinc)
% 0: Dirichlet
% kinc [kx ky] incident wave vector
% This argument is only needed for periodic boundaries.


    if kinc == false
        kinc = [ 0 0 ];
    end
    
    if NGRID(1) == 1
        I = eye();
        DEX = 1i*kinc(1)*I;
        
    else
        n = NGRID(1)*NGRID(2);
        
        if BC(1) == 0
            DEX = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            diagonal_higher(1:NGRID(1):NGRID(1)*NGRID(2)) = 0;
            DEX = spdiags([-diagonal_0th diagonal_higher], [0 1], DEX);
            DEX = DEX/RES(1);
        elseif BC(1) == -2
            DEX = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            diagonal_higher(1:NGRID(1):NGRID(1)*NGRID(2)) = 0;
            diagonal_lower = zeros(n,1);
            periodicityX = NGRID(1)*RES(1);
            diagonal_lower(1:NGRID(1):NGRID(1)*NGRID(2)) = exp(1i*kinc(1)*periodicityX);  
            DEX = spdiags([diagonal_lower -diagonal_0th diagonal_higher], [-NGRID(1)+1 0 1], DEX);
            DEX = DEX/RES(1);
        end  
        
    end

    if NGRID(2) == 1
        I = eye();
        DEY = 1i*kinc(2)*I;

    else
       n = NGRID(1)*NGRID(2);
       
       if BC(2) == 0
            DEY = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            DEY = spdiags([-diagonal_0th diagonal_higher], [0 NGRID(1)], DEY);
            DEY = DEY/RES(2);
        elseif BC(2) == -2
            DEY = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            diagonal_lower = zeros(n,1);
            periodicityY = NGRID(2)*RES(2);
            diagonal_lower(1:1:NGRID(1)*NGRID(2)) = exp(1i*kinc(2)*periodicityY);  
            DEY = spdiags([diagonal_lower -diagonal_0th diagonal_higher], [-(NGRID(1)*NGRID(2)-(NGRID(1)-1))+1 0 NGRID(1)], DEY);
            DEY = DEY/RES(2);
       end  
       
    end
    
    DHX = -(DEX)';
    DHY = -(DEY)';
    
end

function [sx,sy] = calcpml2d(NGRID,NPML)
% CALCPML2D Calculate the PML parameters on a 2D grid
%
% [sx,sy] = calcpml2d(NGRID,NPML);
%
% This MATLAB function calculates the PML parameters sx and sy
% to absorb outgoing waves on a 2D grid.
%
% Input Arguments
% =================
% NGRID Array containing the number of points in the grid
% = [ Nx Ny ]
% NPML Array containing the size of the PML at each boundary
% = [ Nxlo Nxhi Nylo Nyhi ]
%
% Output Arguments
% =================
% sx,sy 2D arrays containing the PML parameters on a 2D grid
    
    a_max = 3;
    p = 3;
    sigmaprime_max = 1;
    eta0 = 376.73032165;

    sx = ones(NGRID);
    sy = sx ;
    
    for nx = 1:NPML(1)
        sx(NPML(1)-nx+1,:) = (1 + a_max*(nx/NPML(1))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*nx)/(2*NPML(1))))^2));

    end

    for nx = 1 : NPML(2)
        sx(NGRID(1)-(2)+nx,:) = (1 + a_max*(nx/NPML(2))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*nx)/(2*NPML(2))))^2));

    end

    for ny = 1:NPML(3)
        sy(:,NPML(3)-ny+1) = (1 + a_max*(ny/NPML(3))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*ny)/(2*NPML(3))))^2));

    end

    for ny = 1 : NPML(4)
        sy(:,NGRID(2)-NPML(4)+ny) = (1 + a_max*(ny/NPML(4))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*ny)/(2*NPML(4))))^2));

    end
end


