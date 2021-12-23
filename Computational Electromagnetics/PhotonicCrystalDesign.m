clc;
clear;

micrometers = 1;
nanometers = micrometers/1000;

LAMBDA = 100*nanometers;

a = 0.125*micrometers; %lattice constant
r = 0.35*a; %hole radius
NP1 = 2; %lattice periods before PC
NP2 = 10; %lattice periods of PC
NP3 = 2; %lattice periods after of PC
NPx = 10;

ur_1 = 1; %relative permeability of grating
er_1 = 1; %dielectric constant of grating
ur_2 = 1.0;
er_2 = 3.5^2;

nL = sqrt(ur_1*er_1);
nH = sqrt(ur_2*er_2);

nmax = max([nH nL]);
NRES = 20;
NPML = 20;
bufy1 = 1*micrometers;
bufy2 = 1*micrometers;

dx = min([LAMBDA])/nmax/NRES;
dy = dx;

nx = ceil(a/dy); dx = a/nx;
ny = ceil(a/dy); dy = a/ny;

NP = NP1 + NP2 + NP3;
Sx = a;
Sy = bufy1 + NP*a + bufy2;
Nx = round(Sx/dx);
Ny = round(Sy/dy) + 2*NPML;

Nx = 2*round(Nx/2)+1;
dx = Sx/Nx;


Nx2 = 2*Nx;
Ny2 = 2*Ny;

dx2 = dx/2;
dy2 = dy/2;

xa = [0:Nx-1]*dx;
xa = xa - mean(xa);
ya = [0:Ny-1]*dy;
xa2 = [0:Nx2-1]*dx2;
xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2;


% Unit Cell
Ny2_uc = round(a/dy2);
UC = zeros(Nx2,Ny2_uc);
ya2_uc = [0:Ny2_uc-1]*dy2;
ya2_uc = ya2_uc - mean(ya2_uc);

[Y,X] = meshgrid(ya2_uc, xa2);

UC = (X.^2 + Y.^2) <= r^2;
UC = nH*(UC == 0) + nL*(UC == 1);

% Photonic Crystal

N2X = nL*ones(Nx2, Ny2);
ny1 = 2*NPML + round(bufy1/dy2);
for np = NP1+1:NP1+NP2
    nya = ny1 + (np - 1)*Ny2_uc;
    nyb = nya + Ny2_uc - 1;
    N2X(:,nya:nyb) = UC;
end
N2X(:,nyb+1:Ny2) = nL;

% Show Device

xa3 = [0:Nx2*NPx-1]*dx2;
P = zeros(Nx2*NPx, Ny2); 
for npx = 1:NPx
    nx1 = (npx-1)*Nx2 + 1;
    nx2 = nx1 + Nx2 - 1;
    P(nx1:nx2, :) = N2X;
end

h = imagesc(xa3/micrometers,ya2/micrometers, P.'); %shows the cartesian plane in the graph but not the ticks
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\mu_r');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;