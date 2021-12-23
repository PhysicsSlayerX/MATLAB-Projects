
% BASED from the University of Texas - El Paso EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This code was written by Francis S. Dela Cruz, and was heavily referenced
% from the lectures on CEM (Computational Electromagnetics) by Dr. Raymond Rumpf. 
% The code was built from scratch with just the help of CEM Lectures.
% The basis of this code is from the Maxwell's Equation that was translated
% into matrices to be solved in MATLAB. This is a benchmarking code for
% future FDTD1D Calculations. Note that this was written with just a basic
% knowledge in MATLAB. Optimization and refactoring of codes are necessary
% to keep the runtime lower and simulation speed faster, as the coder
% implements better code syntax and alogrithms.

% Benchmarked as of August 5, 2020

% This MATLAB Program implements the Finite Difference Time Domain Method in 2D (FDTD-2D).

close all;
clc;
clear all;


t1 = clock;
opengl('software');
centimeters = 1;
millimeters = 0.1*centimeters;
meters = 100*centimeters;
seconds = 1;
minutes = 60*seconds;
c0 = 299792458*meters/seconds;
hertz = 1/seconds;
gigahertz = 1e9*hertz;
e0 = 8.854*10^-12;
%% Simulation Parameters

freq = 5*gigahertz;
fmax = max(freq);
lambda_0 = c0/freq;
er_r = 1;
ur_r = 1;
er_src = 1;
ur_src = 1;
er_t = 1;
ur_t = 1;
ur_bc = 1;
er_bc = 1;
% Material Parameters
ur = 1.0;
er = 3.0;         
ur_rl = 1;
er_rl = 1;
sigma_x = 0; %Testing Conductivity values.
sigma_y = 0; %Testing Conductivity values.
sigma_z = 0; %Testing Conductivity values.

L = 100;
h = 100;
radius = 20;
th = 20*centimeters;
d = 20*centimeters;

NRES = [10 4];
BUFZ = 0.5*[lambda_0 lambda_0 lambda_0 lambda_0]; %Spacer Region
NPML = [20 20 20 20]; %[x-lo x-hi y-lo y-hi]

numgraph_row = 2;
numgraph_col = 2;

freq_fine = 1; % make the frequency resolution finer
counter = 1; % Timestep controller - just for tests

% Movie Creation
MAKE_MOVIE = 0;
movie_title = 'FDTD2D.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    vidObj.FrameRate = 30;
    VidObj.Quality = 100;
    open(vidObj);
end

%% Consideration in Wavelength and Mechanical Resolution

% Wavelength
nref = sqrt(ur_r*er_r);
ntrn = sqrt(ur_t*er_t);
ngrid = sqrt(ur_bc*er_bc);
ndev = sqrt(ur*er);
nreflayer = sqrt(ur_rl*er_rl);
nsrc = sqrt(ur_src*er_src);

nset = [nref ntrn ngrid ndev nsrc nreflayer];

nmax = max(nset);
lambda_min = c0/(fmax*nmax);
delta_lambda = lambda_min/NRES(1);

% Mechanical Features
mech_featX = [L];
mech_featY = [h];
mech_featZ = [th];

dminX = min(mech_featX);
dminY = min(mech_featY);
dminZ = min(mech_featZ);
dmin = min([dminX dminY dminZ]);
delta_d = dmin/NRES(2);

delta_x = min([delta_lambda delta_d]);
delta_y = delta_x;
delta_z = delta_x;

Mx = ceil(dminX/delta_x);
My = ceil(dminY/delta_y);
Mz = ceil(dminZ/delta_z);

dx = dminX/Mx;
dy = dminY/My;
dz = dminZ/Mz;

Nx = round(L/dz) + ceil((2*BUFZ(1) + BUFZ(2))/(nref*dx)) + 3; 
Ny = round(h/dz) + ceil((2*BUFZ(3) + BUFZ(4))/(nref*dy)) + 3; 
%Nz = round(th/dz) + ceil((2*BUFZ(1) + BUFZ(2))/(nref*dz)) + 3; 


% 2X Grid
Nx2 = 2*Nx;
Ny2 = 2*Ny;
dx2 = dx/2;
dy2 = dy/2;

%% Initialize Materials to Free Space

xa2 = [0:Nx2-1]*dx2;
ya2 = [0:Ny2-1]*dy2;
xa2 = xa2 - mean(xa2);
ya2 = ya2 - mean(ya2);

[Y2, X2] = meshgrid(ya2, xa2);

ER2_ = ones(Nx2, Ny2);
UR2_ = ones(Nx2, Ny2);

ER = (X2.^2 + Y2.^2) <= radius^2;
%UR = (X2.^2 + Y2.^2) <= radius^2;
UR = UR2_;

ER = er.*ER + not(ER);
% UR = ur.*UR + not(UR);
%Smoothing Geometry

SmoothFunc = exp(-(X2.^2 + Y2.^2)/0.05^2);
ER2 = fft2(ER).*fft2(SmoothFunc)/sum(SmoothFunc(:));
ER2 = ifftshift(real(ifft2(ER2)));
%ER = ER > 1;
% 
% UR2 = fft2(UR).*fft2(SmoothFunc)/sum(SmoothFunc(:));
% UR2 = ifftshift(real(ifft2(UR2)));
UR2 = UR;
%Extract Grid Parameters


ERxx = ER2(2:2:Nx2, 1:2:Ny2);
ERyy = ER2(1:2:Nx2, 2:2:Ny2);
ERzz = ER2(1:2:Nx2, 1:2:Ny2);

URxx = UR2(1:2:Nx2, 2:2:Ny2);
URyy = UR2(2:2:Nx2, 1:2:Ny2);
URzz = UR2(2:2:Nx2, 2:2:Ny2);

%Plot ot the 2X Grid


fig1 = figure('Color','w');
set(fig1, 'Position', [1 41 1920 963]);
set(fig1, 'Name', '2X Grid');
set(fig1, 'NumberTitle', 'off');

subplot(numgraph_row,numgraph_col,1);

h = imagesc(xa2,ya2, ER2.');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\epsilon_r');
axis([min(xa2) max(xa2) min(ya2) max(ya2)]);
axis equal tight;
colorbar;

subplot(numgraph_row,numgraph_col,3);

g = imagesc(xa2,ya2, UR2.',[1 10]);
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\mu_r');
axis([min(xa2) max(xa2) min(ya2) max(ya2)]);
axis equal tight;
colorbar;

drawnow; 

% Compute Time Step
delta_min = min([dx dy dz]);
dt = (ngrid*delta_min)/(2*c0);

%% Compensate for Numerical Dispersion

% k0 = (2*pi*fmax)./c0;
% n_avg = mean(nset);
% numdisp = c0*dt/(n_avg*dz) *sin(k0*n_avg*dz/2)/sin(c0*k0*dt/2);
% 
% UR = numdisp*UR;
% ER = numdisp*ER;


%% Reflection, Source, and Transmission Layers


er_ref = ER2(:, 2*NPML(3) + 10); %One layer in the reflection region's relative permittivity 
er_ref = mean(er_ref(:)); %mean relative permittivity of the reflection grid layer
er_trn = ER2(:, Ny2 - 2*NPML(3) - 5);%One layer in the transmission region's relative permittivity
er_trn = mean(er_trn(:)); %mean relative permittivity of the transmission grid layer
er_src = ER2(:, 2*NPML(3) + 11); %One layer AFTER the reflection region's relative permittivity
er_src = mean(er_src(:)); %mean relative permittivity of the source grid layer

ur_ref = UR2(:, 2*NPML(4) + 10); %One layer in the reflection region's relative permeability 
ur_ref = mean(ur_ref(:)); %mean relative permittivity of the grid layer
ur_trn = UR2(:, Ny2 - 2*NPML(4) - 5);%One layer in the transmission region's relative permeability
ur_trn = mean(ur_trn(:)); %mean relative permittivity of the grid layer
ur_src = UR2(:, 2*NPML(3) + 11); %One layer AFTER the reflection region's relative permittivity
ur_src = mean(ur_src(:)); %mean relative permittivity of the source grid layer

nref = sqrt(er_ref*ur_ref); %refractive index of the reflection region
ntrn = sqrt(er_trn*ur_trn); %refractive index of the transmission region
nsrc = sqrt(er_src*ur_src); %refractive index of the source region

if er_ref <0 && ur_ref<0
    nref = - nref;
end

if er_trn <0 && ur_trn <0
    ntrn = - ntrn;
end



%% Compute Source

tau = 0.5/fmax;
time0 = 6*tau;

Tprop = (nmax*Ny*dy)/c0;
Ttotal = 12*tau + 5*Tprop;
STEPS = ceil(Ttotal/dt);

Time_element = [0:STEPS-1]*dt;
% Ez Mode
nysrc =  2*NPML(3) + 11;
Ezsrc = exp(-((Time_element - time0)/tau).^2);
A = sqrt(er_src/ur_src); %-((sigma_yy/(fmax*e0*sqrt(ur_src*er_src)))+sqrt(er_src/ur_src)); -> This is a loss term that needs to be edited
sigma_t = (nsrc*dy)/(2*c0)+dt/2;
Hxsrc= A.*exp(-((Time_element - time0 + sigma_t)/tau).^2);

% fig2 = figure('Color','w');
% 
% subplot(1, 2, 1);
% set(fig2, 'Position', [1 41 1920 963]);
% set(fig2, 'Name', 'FDTD-2D Analysis - Source Function');
% set(fig2, 'NumberTitle', 'off');
% x_axis = Time_element;
% plot(x_axis, Ezsrc,'-b','LineWidth', 3); hold on;
% plot(x_axis, Hxsrc,'-r','LineWidth', 3); hold on;
% axis([min(Time_element) max(Time_element) -1 1]);
% xlabel('Time (s)');
% ylabel('Amplitude','Rotation',90);
% title('Source Functions');
% h = legend('Ey','Hx','Location','NorthEastOutside');
% set(h,'LineWidth',2);


%Hz Mode
Exsrc = exp(-((Time_element - time0)/tau).^2);
A = -sqrt(er_src/ur_src); %-((sigma_yy/(fmax*e0*sqrt(ur_src*er_src)))+sqrt(er_src/ur_src)); -> This is a loss term that needs to be edited
sigma_t = (nsrc*dy)/(2*c0)+dt/2;
Hzsrc = A.*exp(-((Time_element - time0 + sigma_t)/tau).^2);

% 
% subplot(1,2,2);
% plot(x_axis, Exsrc,'-b','LineWidth', 3); hold on;
% plot(x_axis, Hzsrc,'-r','LineWidth', 3); hold off;
% axis([min(Time_element) max(Time_element) -1 1]);
% xlabel('Time (s)');
% ylabel('Amplitude','Rotation',90);
% title('Source Functions');
% h = legend('Ey','Hx','Location','NorthEastOutside');
% set(h,'LineWidth',2);



% Initialize Fourier Transform
% if freq_fine == 1
%     NFREQ = STEPS;
% else
%     NFREQ = 100;
% end
% FREQ = linspace(0, freq, NFREQ);
% K = exp(-1i*2*pi*dt.*FREQ);
% REFe = zeros(1, NFREQ);
% TRNe = zeros(1, NFREQ);
% REFh = zeros(1, NFREQ);
% TRNh = zeros(1, NFREQ);
% SRC = zeros(1, NFREQ);

%% Compute for Perfectly Matched Layer Parameters

sigx = zeros(Nx2, Ny2);
for nx = 1: 2*NPML(1)
    nx1 = 2*NPML(1) - nx +1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1: 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end

sigy = zeros(Nx2, Ny2);
for ny = 1: 2*NPML(3)
    ny1 = 2*NPML(3) - ny +1;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1: 2*NPML(2)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end


subplot(numgraph_row,numgraph_col,2);

img = imagesc(xa2,ya2, sigx.');
img2 = get(img, 'Parent');
set(img2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\sigma_x');
axis([min(xa2) max(xa2) min(ya2) max(ya2)]);
axis equal tight;
colorbar;

subplot(numgraph_row,numgraph_col,4);

img_2 = imagesc(xa2,ya2, sigy.');
img_3 = get(img_2, 'Parent');
set(img_3,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\sigma_y');
axis([min(xa2) max(xa2) min(ya2) max(ya2)]);
axis equal tight;
colorbar;

drawnow; 

%% Compute for Update Coefficients

%Ez Mode
sigHx = sigx(1:2:Nx2, 2:2:Ny2);
sigHy = sigy(1:2:Nx2, 2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = - c0./URxx./mHx0;
mHx3 = -(c0*dt/e0)*sigHx./URxx./mHx0;

sigHx = sigx(2:2:Nx2, 1:2:Ny2);
sigHy = sigy(2:2:Nx2, 1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./URyy./mHy0;
mHy3 = -(c0*dt/e0)*sigHy./URyy./mHy0;

sigDx = sigx(1:2:Nx2, 1:2:Ny2);
sigDy = sigy(1:2:Nx2, 1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1./mDz0;
mDz2 = c0./mDz0;
mDz4 = -(dt/e0^2)*sigDx.*sigDy./mDz0;

mEz1 = 1./ERzz;

%Hz Mode

sigDx = sigx(2:2:Nx2, 1:2:Ny2);
sigDy = sigy(2:2:Nx2, 1:2:Ny2);
mDx0 = (1/dt) + sigDy/(2*e0);
mDx1 = ((1/dt) - sigDy/(2*e0))./mDx0;
mDx2 = c0./mDx0;
mDx3 = (c0*dt/e0)*sigDx./mDx0;

sigDx = sigx(1:2:Nx2, 2:2:Ny2);
sigDy = sigy(1:2:Nx2, 2:2:Ny2);
mDy0 = (1/dt) + sigDx/(2*e0);
mDy1 = ((1/dt) - sigDx/(2*e0))./mDy0;
mDy2 = c0./mDy0;
mDy3 = (c0*dt/e0)*sigDy./mDy0;

sigHx = sigx(2:2:Nx2, 2:2:Ny2);
sigHy = sigy(2:2:Nx2, 2:2:Ny2);
mHz0 = (1/dt) + (sigHx + sigHy)/(2*e0) + sigHx.*sigHy*(dt/4/e0^2);
mHz1 = (1/dt) - (sigHx + sigHy)/(2*e0) - sigHx.*sigHy*(dt/4/e0^2);
mHz1 = mHz1./mHz0;
mHz2 = - c0./URzz./mHz0;
mHz4 = -(dt/e0^2)*sigHx.*sigHy./mHz0;

mEx1 = 1./ERxx;
mEy1 = 1./ERyy;


%% Field Initialization

%Ez Mode

Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);
CEx = zeros(Nx, Ny);
CEy = zeros(Nx, Ny);
CHx = zeros(Nx, Ny);
ICEx = zeros(Nx, Ny);
ICEy = zeros(Nx, Ny);
IDz = zeros(Nx, Ny);

%Hz mode

Ex = zeros(Nx, Ny);
Ey = zeros(Nx, Ny);
Dx = zeros(Nx, Ny);
Dy = zeros(Nx, Ny);
Hz = zeros(Nx, Ny);
CHy = zeros(Nx, Ny);
CHz = zeros(Nx, Ny);
CEz = zeros(Nx, Ny);
ICHx = zeros(Nx, Ny);
ICHy = zeros(Nx, Ny);
IHz = zeros(Nx, Ny);

%% Main FDTD Loop
fig2 = figure('Color','w');
set(fig2, 'Position', [1 41 1920 963]);
set(fig2, 'Name', 'FDTD-2D');
set(fig2, 'NumberTitle', 'off');

for T = 1:counter:STEPS
%Ez Mode
    %Compute for CEx
    % Dirichlet Boundary Condition for CEx

    for nx = 1:Nx
       for ny = 1 : Ny - 1
          CEx(nx, ny) = (Ez(nx, ny+1) - Ez(nx,ny))/dy; 
       end
       CEx(nx, Ny) = (0 - Ez(nx,ny))/dy;
    end
    %Compute for CEy
    % Dirichlet Boundary Condition for CEy

    for ny = 1:Ny
       for nx = 1 : Nx - 1
          CEy(nx, ny) = -(Ez(nx+1, ny) - Ez(nx,ny))/dx; 
       end
       CEy(Nx, ny) = - (0 - Ez(Nx,ny))/dx;
    end
    %Inject TF/SF Source into Curl of E
    for nx = 1 : Nx
        CEx(nx, nysrc) = CEx(nx, nysrc) - Ezsrc(T)/dy;
    end
    %Update H Integrations
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    %Update H Field
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
    
    % Dirichlet Boundary Condition for CHz

    CHz(1,1) = (Hy(1,1) - 0)/dx...
             - (Hx(1,1) - 0)/dy;
    for nx = 2 : Nx
        CHz(nx,1) = (Hy(nx,1) - Hy(nx-1, 1))/dx...
                  - (Hx(nx,1) - 0)/dy;
    end
    for ny = 2 : Ny
        CHz(1,ny) = (Hy(1,ny) - 0)/dx...
                  - (Hx(1,ny) - Hx(1, ny-1))/dy;
        for nx = 2 : Nx
            CHz(nx,ny) = (Hy(nx,ny) - Hy(nx-1, ny))/dx...
                      - (Hx(nx,ny) - Hx(nx, ny-1))/dy;
        end
    end
    %Inject TF/SF Source into Curl of E
    for nx = 1 : Nx
        CHz(nx, nysrc) = CHz(nx, nysrc) + Hxsrc(T)./dy;
    end
    % Update D Integrations
    IDz = IDz + Dz;
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
   
    %Update Ez
    Ez = mEz1.*Dz;
%Hz Mode
    %Dirichlet Boundary Condition for CEz
    CEz(Nx, Ny) = (0 - Ey(Nx, Ny))/dx - (0 - Ex(Nx, Ny))/dy;
    for nx = 1: Nx-1
       CEz(nx, Ny) = (Ey(nx+1, Ny) - Ey(nx, Ny))/dx - (0 - Ex(nx, Ny))/dy; 
    end
    for ny = 1: Ny-1
       CEz(Nx, ny) = (0 - Ey(Nx, ny))/dx - (Ex(Nx, ny+1) - Ex(Nx, ny))/dy;
       for nx = 1: Nx-1
           CEz(nx, ny) = (Ey(nx+1, ny) - Ey(nx, ny))/dx - (Ex(nx, ny+1) - Ex(nx, ny))/dy;
       end
    end
    %Inject TF/SF Source into Curl of E
    for nx = 1 : Nx
        CEz(nx, nysrc) = CEz(nx, nysrc) + Exsrc(T)./dy;
    end
    %Update H Integrations
    IHz = IHz + Hz;
    % Update H-Field
    Hz = mHz1.*Hz + mHz2.*CEz + mHz4.*IHz;

    %Dirichlet Boundary Condition for CHx
    for nx = 1: Nx
        CHx(nx, 1) = (Hz(nx, 1) - 0)/dy;
        for ny = 2 : Ny
          CHx(nx, ny) = (Hz(nx, ny) - Hz(nx, ny-1))/dy; 
        end
    end    
    
    %Dirichlet Boundary Condition for CHy
    
    for ny = 1: Ny
       CHy(1, ny) = -(Hz(1, ny) - 0)/dx;
       for nx = 2 : Nx
           CHy(nx, ny) = -(Hz(nx, ny) - Hz(nx - 1, ny))/dx;
       end
    end
    %Inject TF/SF Source into Curl of E
    for nx = 1 : Nx
        CHx(nx, nysrc) = CHx(nx, nysrc) - Hzsrc(T)./dy;
    end
    % Update D Integrations
    
    ICHx = ICHx + CHx;
    ICHy = ICHy + CHy;
    
    %Update D Field
    
    Dx = mDx1.*Dx + mDx2.*CHx + mDx3.*ICHx;
    Dy = mDy1.*Dy + mDy2.*CHy + mDy3.*ICHy;
    
    % Update E Field
    
    Ex = mEx1.*Dx;
    Ey = mEy1.*Dy;
    
    %Update Fourier Transforms
%     for nf = 1: NFREQ
%         REFe(nf) = REFe(nf) + (K(nf)^T)*Ey(nzref);
%         TRNe(nf) = TRNe(nf) + (K(nf)^T)*Ey(nztrn);
%         REFh(nf) = REFh(nf) + (K(nf)^T)*Hx(nzref);
%         TRNh(nf) = TRNh(nf) + (K(nf)^T)*Hx(nztrn);
%         SRC(nf) = SRC(nf) + (K(nf)^T)*Eysrc(T);
%     end
    if MAKE_MOVIE == 1
        clf;
    end
% Ez Mode    
    % Update Visualization
    subplot(1, 2, 1);
    x_axis = [0:Nx-1]*dx;
    y_axis = [0:Ny-1]*dy;
    x_axis = x_axis - mean(x_axis);
    y_axis = y_axis - mean(y_axis); 
    imagesc(x_axis, y_axis, Ez.');
    axis([min(x_axis) max(x_axis) min(y_axis) max(y_axis) ]);
    xlabel('x (cm)');
    ylabel('y (cm)','Rotation',90);
    title(['Timestep ' num2str(T) ' of ' num2str(ceil(Ttotal/dt)-1)], 'Interpreter', 'LaTex', 'FontSize', 18);
%     h = legend('Anti-reflection Layer', 'Device', 'Anti-reflection Layer','Ey','Hx','Location','NorthEastOutside');
%     set(h,'LineWidth',2);
    colorbar;
    caxis([0 1]);
    axis equal tight;
    subplot(1, 2, 2);
    
    imagesc(x_axis, y_axis, abs(Hz).');
    axis([min(x_axis) max(x_axis) min(y_axis) max(y_axis) ]);
    xlabel('x (cm)');
    ylabel('y (cm)','Rotation',90);
    title(['Timestep ' num2str(T) ' of ' num2str(ceil(Ttotal/dt)-1)]);
%     h = legend('Anti-reflection Layer', 'Device', 'Anti-reflection Layer','Ey','Hx','Location','NorthEastOutside');
%     set(h,'LineWidth',2);
    colorbar;
    axis equal tight;
    caxis([0 1]);
    
    
%     x_axis = FREQ;
%     plot(x_axis, (abs(REFe./SRC).^2),'-r','LineWidth', 3); hold on;
%     plot(x_axis, (abs(TRNe./SRC).^2),'-b', 'LineWidth', 3);
%     plot(x_axis, (abs(REFe./SRC).^2)+(abs(TRNe./SRC).^2),'--k', 'LineWidth', 3);
%     plot(x_axis, (1 - (abs(REFe./SRC).^2)+(abs(TRNe./SRC).^2)),'--g', 'LineWidth', 3); hold off;
%     axis([min(FREQ) max(FREQ) 0 1.25]);
%     xlabel('Frequency (GHz)');
%     ylabel('Reflectance and Transmittance','Rotation',90);
%     title('Spectral Response - Ey');
%     h = legend('REFe','TRNe','CONe','ABSe','Location','NorthEastOutside');
%     set(h,'LineWidth',2);
    drawnow;
    % Add Frame to AVI
    if MAKE_MOVIE == 1
        Frames = getframe(fig2);
        writeVideo(vidObj, Frames);
    end
end

if MAKE_MOVIE == 1
        close(vidObj);
end
% Finish Fourier Transform
% REFe = REFe*dt;
% TRNe = TRNe*dt;
% REFh = REFh*dt;
% TRNh = TRNh*dt;
% SRC = SRC*dt;
% 
% % Compute Reflectance and Transmittance
% 
% REFe = abs(REFe./SRC).^2;
% TRNe = abs(TRNe./SRC).^2;
% CONe = REFe + TRNe;
% 
% REFh = abs(REFh./SRC).^2;
% TRNh = abs(TRNh./SRC).^2;
% CONh = REFh + TRNh;
% 
% ABSe = 1 - CONe;
% ABSh = 1 - CONh;
% 
% RdBe = 10*log10(REFe);
% TdBe = 10*log10(TRNe);
% CdBe = 10*log10(CONe);
% 
% RdBh = 10*log10(REFh);
% TdBh = 10*log10(TRNh);
% CdBh = 10*log10(CONh);


% Visualization
% fig2 = figure('Color','w');
% set(fig2, 'Position', [1 41 1920 963]);
% set(fig2, 'Name', 'FDTD-1D Analysis - Source Function');
% set(fig2, 'NumberTitle', 'off');
% x_axis = Time_element;
% plot(x_axis, Eysrc,'-b','LineWidth', 3); hold on;
% plot(x_axis, Hxsrc,'-r','LineWidth', 3); hold off;
% axis([min(Time_element) max(Time_element) -1 1]);
% xlabel('Time (s)');
% ylabel('Amplitude','Rotation',90);
% title('Source Functions');
% h = legend('Ey','Hx','Location','NorthEastOutside');
% set(h,'LineWidth',2);
%  
% fig3 = figure('Color','w');
% set(fig3, 'Position', [1 41 1920 963]);
% set(fig3, 'Name', 'FDTD-1D Analysis - Spectral Response');
% set(fig3, 'NumberTitle', 'off');
% 
% subplot(2, 2, 1);
% 
% x_axis = FREQ;
% plot(x_axis, REFe,'-r', 'LineWidth', 2.5); hold on;
% plot(x_axis, TRNe,'-b','LineWidth', 2.5);
% plot(x_axis, CONe,'--k','LineWidth', 2.5); 
% plot(x_axis, ABSe,'--g','LineWidth', 2.5); hold off;
% axis([min(FREQ) max(FREQ) 0 1.25]);
% xlabel('Frequency (GHz)');
% ylabel('Reflectance and Transmittance','Rotation',90);
% title('Spectral Response - Ey');
% h = legend('REFe','TRNe','CONe','Location','NorthEastOutside');
% set(h,'LineWidth',2);
% 
% subplot(2, 2, 2);
% 
% x_axis = FREQ;
% plot(x_axis, RdBe,'-r','LineWidth', 2.5); hold on;
% plot(x_axis, TdBe,'-b','LineWidth', 2.5);
% plot(x_axis, CdBe,'--k','LineWidth', 2.5); hold off;
% axis([min(FREQ) max(FREQ) -30 0]);
% xlabel('Frequency (GHz)');
% ylabel('Reflectance and Transmittance (dB)','Rotation',90);
% title('Spectral Response - Ey(dB)');
% h = legend('REFe(dB)','TRNe(dB)','CONe(dB)','Location','NorthEastOutside');
% set(h,'LineWidth',2);
% 
% subplot(2, 2, 3);
% 
% x_axis = FREQ;
% plot(x_axis, REFh,'-r','LineWidth', 2.5); hold on;
% plot(x_axis, TRNh,'-b','LineWidth', 2.5);
% plot(x_axis, CONh,'--k','LineWidth', 2.5); 
% plot(x_axis, ABSh,'--g','LineWidth', 2.5); hold off;
% axis([min(FREQ) max(FREQ) 0 1.25]);
% xlabel('Frequency (GHz)');
% ylabel('Reflectance and Transmittance','Rotation',90);
% title('Spectral Response - Hx');
% h = legend('REFh','TRNh','CONh','Location','NorthEastOutside');
% set(h,'LineWidth',2);
% 
% subplot(2, 2, 4);
% 
% x_axis = FREQ;
% plot(x_axis, RdBh,'-r','LineWidth', 2.5); hold on;
% plot(x_axis, TdBh,'-b','LineWidth', 2.5);
% plot(x_axis, CdBh,'--k','LineWidth', 2.5); hold off;
% axis([min(FREQ) max(FREQ) -30 0]);
% xlabel('Frequency (GHz)');
% ylabel('Reflectance and Transmittance (dB)','Rotation',90);
% title('Spectral Response - Hx(dB)');
% h = legend('REFh(dB)','TRNh(dB)','CONh(dB)','Location','NorthEastOutside');
% set(h,'LineWidth',2);

t2 = clock;
t = etime(t2,t1);
disp(['Elapsed time is ' num2str(t/minutes) ...
    ' minutes.']);
