% Homework #7, Problem 1
% EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This MATLAB script file implements the FDFD method
% to model transmission and reflection from a grating.
% INITIALIZE MATLAB
close all;
clc;
clear all;

% OPEN FIGURE WINDOW
fig = figure('Color','w');
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
t = 0.1510*lamd; %substrate thickness
mech_feat = [x1 x2 x3 L d t];
ur = 1.0; %relative permeability of grating
er = 9.0; %dielectric constant of grating
% EXTERNAL MATERIALS
ur1 = 1.0; %permeability in the reflection region
er1 = 1.0; %permittivity in the reflection region
ur2 = 1.0; %permeability in the transmission region
er2 = 9.0; %permittivity in the transmission region
% GRID PARAMETERS
NRES = 40; %grid resolution
BUFZ = 2*lam0 * [1 1]; %spacer region above and below grating
NPML = [0 0 20 20]; %size of PML at top and bottom of grid
xbc = -2;
ybc = 0;

% OPTIMIZED PARAMETERS

% Consider the Wavelength
n_device = sqrt(ur*er);
nref = sqrt(ur1*er1);
ntrn = sqrt(ur2*er2);

n_max = max([n_device nref ntrn]);

lambda_min = lam0/n_max;

delta_gamma = lambda_min/NRES;

% Consider the Mechanical Features

d_min = min(mech_feat);

delta_d = d_min/6;

% Choosing which is smaller between the parameters
delta_x = min(delta_gamma, delta_d);
delta_y = delta_x;

% Resolving Critical Dimensions through "Snapping"
xTdev = x1+x2+x3;
yTdev = d+t;

Mx = ceil(xTdev/delta_x);
My = ceil(yTdev/delta_y);

% Total Grid Size

dx = xTdev/Mx;
dy = yTdev/My;

Nx = round(L/dx);
Nx = 2*round(Nx/2)+1;
Ny = round(yTdev/dy)+ 2.*NPML(3)+ 2.*ceil(2*lam0/(nref*dy));
%

%The 2x Grid Parameters

Nx2 = 2*Nx;
Ny2 = 2*Ny;

dx2 = dx/2;
dy2 = dy/2;

% Building the Simulation Matrix and the device

nx1 = 1+floor((Nx2 - round(xTdev/dx2))/2);
nx2 = nx1 + round(x1/dx2) - 1;
nx3 = nx2+1;
nx4 = nx3 + round(x2/dx2) - 1;
nx5 = nx4+1;
nx6 = nx5 + round(x3/dx2) - 1;

ny1 = 1+floor((Ny2 - round(yTdev/dy2))/2);
ny2 = ny1 + round(d/dy2) - 1;
ny3 = ny2+1;
ny4 = ny3 + round(t/dy2) - 1;

ER2 = er1*ones(Nx2,Ny2);
ER2(:, ny1:ny2) = er;
%ER2(nx5:nx6, ny1:ny2) = er;
ER2(:,ny3:ny4) = er;
ER2(:,ny4+1:Ny2) = er2;

UR2 = ur1*ones(Nx2,Ny2);

% ER2(:, 2*NPML(3) + 2)= 5;
% ER2(:, 2*NPML(3) + 5)= 6;
% ER2(:, Ny2 - 2*NPML(3) - 2) = 7;
% 
% UR2(:, 2*NPML(4) + 2) = 5;
% UR2(:, 2*NPML(4) + 5) = 6;
% UR2(:, Ny2 - 2*NPML(4) - 2) = 7;

% The Finite Difference Frequency Domain Method

% Material Properties in the Reflected and Transmitted Regions
k0 = 2*pi./lam0;
er_ref = ER2(:, 2*NPML(3) + 10);
er_ref = mean(er_ref(:));
er_trn = ER2(:, Ny2 - 2*NPML(3) - 5);
er_trn = mean(er_trn(:));

ur_ref = UR2(:, 2*NPML(4) + 10);
ur_ref = mean(ur_ref(:));
ur_trn = UR2(:, Ny2 - 2*NPML(4) - 5);
ur_trn = mean(ur_trn(:));
nref(er_ref*ur_ref);

NGRID2X = [Nx2 Ny2];
[sx, sy] = calcpml2d(NGRID2X,2*NPML);

% Incorporate the PML to the grid:

URxx = UR2./sx.*sy;
URyy = UR2.*sx./sy;
URzz = UR2.*sx.*sy;
ERxx = ER2./sx.*sy;
ERyy = ER2.*sx./sy;
ERzz = ER2.*sx.*sy;

% Overlay Materials Onto 1x Grids

URxx = URxx(1:2:Nx2,2:2:Ny2);
URyy = URyy(2:2:Nx2,1:2:Ny2);
URzz = URzz(2:2:Nx2,2:2:Ny2);
ERxx = ERxx(2:2:Nx2,1:2:Ny2);
ERyy = ERyy(1:2:Nx2,2:2:Ny2);
ERzz = ERzz(1:2:Nx2,1:2:Ny2);

% Compute the Wave Vector Terms

kinc = k0.*nref.*[sin(theta); cos(theta)];

m = [-floor(Nx/2):floor(Nx/2)]';
kx_m = kinc(1) - m*(2*pi/(Nx*dx));

ky_ref_m = (sqrt((k0*nref)^2 - (kx_m).^2));

ky_trn_m = (sqrt((k0*ntrn)^2 - (kx_m).^2));

% Convert Diagonal Materials Matrices

ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% Construc the Derivative Matrices

NGRID = [Nx Ny];
RES = [dx dy];
BC = [xbc ybc];

[DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);

% Compute the Wave Matrix A
% E-Mode and H-Mode

Ae = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
Ah = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;

% Compute the Source Field
x1a = m*dx;
y1a = [1:Ny]*dy;
[Y1a, X1a] = meshgrid(y1a, x1a);
f_src = exp(1i.*(kinc(1).*X1a + kinc(2).*Y1a));
fsrc = f_src(:);

% Compute the Scattered-Field Masking Matrix, Q

Q = sparse(Nx, Ny);
sr = 1;
Q(:, 1:NPML(3)+20) = sr;
%Q = reshape(Q,[],1);
Q = diag(sparse(Q(:)));

% Compute the source vector, b; 

bE = (Q*Ae - Ae*Q)*fsrc;
bH = (Q*Ah - Ah*Q)*fsrc;

efield = Ae\bE;
efield = full(efield);
efield = reshape(efield,Nx,Ny);

hfield = Ah\bH;
hfield = full(hfield);
hfield = reshape(hfield, Nx, Ny);

% Extract Transmitted and Reflected Fields

Eref = efield(:, NPML(3)+10);
Etrn = efield(:, Ny - NPML(3)-10);

Href = hfield(:, NPML(4) + 10);
Htrn = hfield(:, Ny-NPML(4) -10);

% Remove the Phase Tilt
x_value = [1:Nx]'*dx;
Aref = Eref.*exp(-1i*(kinc(1)*x_value));
Atrn = Etrn.*exp(-1i*(kinc(1)*x_value));

Bref = Href.*exp(-1i*(kinc(1)*x_value));
Btrn = Htrn.*exp(-1i*(kinc(1)*x_value));

% Calculate the Complex Amplitudes of the Spatial Harmonics

Sref = flipud(fftshift(fft(Aref)))/Nx;
Strn = flipud(fftshift(fft(Atrn)))/Nx;
Uref = flipud(fftshift(fft(Bref)))/Nx;
Utrn = flipud(fftshift(fft(Btrn)))/Nx;

% Calculate Diffraction Efficiencies

RE_m = 100.*abs(Sref).^2.*(real(ky_ref_m./ur1)./real(kinc(2)/ur1));
RH_m = 100.*abs(Uref).^2.*(real(ky_ref_m./er1)./real(kinc(2)/er1));

TE_m = 100.*abs(Strn).^2.*(real(ky_trn_m./ur2)./real(kinc(2)/ur1));
TH_m = 100.*abs(Utrn).^2.*(real(ky_trn_m./er2)./real(kinc(2)/er1));


% Reflectance and Transmittance

REF_E = sum(RE_m);
REF_H = sum(RH_m);
TRN_E = sum(TE_m);
TRN_H = sum(TH_m);

CON_E = REF_E + TRN_E;
CON_H = REF_H + TRN_H;

% Printing the Results

disp(['R_E = ' num2str(REF_E)]);
disp(['T_E = ' num2str(TRN_E)]);
disp(['R_H = ' num2str(REF_H)]);
disp(['T_H = ' num2str(TRN_H)]);
disp(['CON_E = ' num2str(CON_E)]);
disp(['CON_H = ' num2str(CON_H)]);


% Plotting the 2D Device
numgraph = 5;

subplot(1,numgraph,2);

xa = [-Nx2/2:Nx2/2]*dx2;
ya = [0:Ny2-1]*dy2;
[Y, X] = meshgrid(ya, xa);

h = imagesc(xa,ya, ER2.');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\epsilon_r');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;
%colormap(jet(1024));

subplot(1,numgraph,1);

g = imagesc(xa,ya, UR2.',[1 10]);
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\mu_r');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
%colormap(jet(1024));

subplot(1,numgraph,3);
xb = [m]*dx;
yb = [1:Ny]*dy;
[Yb, Xb] = meshgrid(yb, xb);
h = imagesc(xb,yb, real(hfield).');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('H-mode at 24 GHz');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;
%colormap(jet(1024));

subplot(1,numgraph,4);

h = imagesc(xb,yb, real(efield).');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('E-mode at 24 GHz');
shading interp;
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;
%colormap(jet(1024));

subplot(1,numgraph,5);

h = imagesc(xb,yb, real(f_src).');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('f_src');
shading interp;
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar('Limits', [-10 10]);
%colormap(jet(1024));


set(gcf, 'Position', [1920 50 1900 870])



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

