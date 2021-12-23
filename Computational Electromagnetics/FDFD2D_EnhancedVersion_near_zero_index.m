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

% START TIMER
t1 = clock;

% OPEN FIGURE WINDOW
opengl('software');
fig = figure('Color','w');
set(gcf, 'Position', [1 41 1920 963]);
set(fig, 'Name', 'FDFD Analysis');
set(fig, 'NumberTitle', 'off');

% FORMAT SIGNIFICANT DIGITS
format short g; 
% UNITS
centimeters = 1;
millimeters = 0.1 * centimeters;
meters = 100 * centimeters;
micrometers = 1;
nanometers = 1/1000;
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
%f0 = [11]* gigahertz; %operating frequency
%LAMBDA = c0./f0; %operating wavelength in free space
LAMBDA = [978]*nanometers;
SRC.theta = [45] * degrees; %angle of incidence
THETA = SRC.theta;
SRC.MODE = 'E'; %electromagnetic mode, 'E' or 'H'
SRC.scalar_pos = 0.30; % 0 < x <= 1 position of the beam.
SRC.sourcetype = 1; % 0,1 for source type. 0-plane, 1-gaussian
SRC.w = min(2*LAMBDA);
% GRATING PARAMETERS
%fd = 8 * gigahertz; %design frequency
%lamd = c0/fd; %design wavelength
%lamd = 1;
L = 10000*nanometers; %length of grating
width = 0.5*L;
height = 500*nanometers;  %tooth height
subsheight = 300*nanometers;
% h = 0.85*LAMBDA(1);  %tooth height
ur = 1; %relative permeability of grating
er = 0.00001^2; %dielectric constant of grating
ur_2 = 1.0;
er_2 = 0.001^2;
ur_3 = 1.0;
er_3 = 10^2;
% EXTERNAL MATERIALS
ur1 = 1.0; %permeability in the reflection region
er1 = 1.0; %permittivity in the reflection region
ur2 = 1.0; %permeability in the transmission region
er2 = 1.0; %permittivity in the transmission region

% GRID PARAMETERS
NRES = [70]; %grid resolution
BUFZ = 2*max(LAMBDA); %spacer region above and below grating
DEV.NPML = [200 200 10 10]; %size of PML at top and bottom of grid
xbc = 0; % xbc tells if the x-axis will have DIrichlet or Pseudo-Periodic/Floquet Boundary Condition
ybc = 0; % xbc tells if the y-axis will have DIrichlet or Pseudo-Periodic/Floquet Boundary Condition
DEV.BC = [xbc ybc];  % Boundary Condition: 0 for Dirichlet, -2 for Pseudo-Periodic/Floquet
%GratingPeriod = 10; % For layered grating simulation

% Plotting Position
numgraph_row = 2;
numgraph_col = 2;
responsegraph = 0;

% Phase Visual
phasegraph = 0;
NPHI = 100;
phase = linspace(0,8*pi, NPHI);

% Movie Creation
MAKE_MOVIE = 0;
movie_title = 'Zero-index_45deg.mp4';
if MAKE_MOVIE == 1
    vidObj = VideoWriter(movie_title, 'MPEG-4');
    open(vidObj);
end

CONV = 0; % Enable Convergence Test
FREQ = 1; % Enable Frequency Sweep
ANGLE = 0; % Enable Angle Sweep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if CONV == 1
    %% Convergence Test
    REF2 = zeros(length(NRES),1);
    TRN2 = zeros(length(NRES),1);
    CON2 = zeros(length(NRES),1);
    for index = 1:length(NRES)
        % Refractive Indices of the Simulation Setup

        n_device = sqrt(ur*er); % refractive index of the device
        n_ref = sqrt(ur1*er1); % refractive index of the reflection region
        n_trn = sqrt(ur2*er2); % refractive index of the transmission region

        % Consider the Wavelength for resolving the resolution of the simulation

        n_max = max([n_device n_ref n_trn]);
        lambda_min = min(LAMBDA)/n_max;
        delta_gammaX = lambda_min/NRES(index);
        delta_gammaY = lambda_min/NRES(index);

        % Consider the Mechanical Features for optimizing the simulaiton
        mech_featX = [width L]; %array of mechanical features
        mech_featY = [height];
        d_minX = min(mech_featX)/4;
        d_minY = min(mech_featY)/4;

        % Choosing which is smaller between the parameters (to resolve smallest
        % features

        delta_x = min(delta_gammaX, d_minX);
        delta_y = min(delta_gammaY, d_minY);

        % Resolving Critical Dimensions through "Snapping"
        xTdev = min(mech_featX); %the critical dimension here is the duty cycle of the grating
        yTdev = min(mech_featY); %the height and the substrate thickness

        Mx = ceil(xTdev/delta_x);
        My = ceil(yTdev/delta_y);

        % Total Grid Size

        dx = xTdev/Mx; %the adjusted dx
        dy = yTdev/My; %the adjusted dy

        Nx = round(L/dx);
        Nx = 2*round(Nx/2)+1; %physical size along x-axis
        Ny = round(yTdev/dy)+ 2.*DEV.NPML(3)+ 2.*ceil(max(2*BUFZ)/(n_ref*dy)); %physical size along y-axis

        %The 2x Grid Parameters

        Nx2 = 2*Nx; %The 2x Grid requires 2x of the number of elements for both x and y
        Ny2 = 2*Ny;

        dx2 = dx/2;
        dy2 = dy/2;

        % Building the Simulation Matrix and the device

        pos_x1 = round((Nx2 - width/dx2)/2);
        pos_x2 = round((Nx2 + width/dx2)/2);

        pos_y1 = (Ny2 - 2*DEV.NPML(1) + 2*round(BUFZ/dy2));
        pos_y1 = 2*round(pos_y1/2)+1;
        pos_y2 = pos_y1 + round(height/dy2) - 1;

        DEV.ER2 = er1*ones(Nx2,Ny2);
        DEV.ER2(:,pos_y1:pos_y2) = er_2;
        DEV.ER2(pos_x1:pos_x2,pos_y1:pos_y2) = er; 
        DEV.ER2(:, pos_y2+1:Ny2) = er2;

        DEV.UR2 = ur1*ones(Nx2,Ny2); %Since the relative permeability remains 1 (non-magnetic), fill the entire grid with ur1 = 1
        DEV.UR2(:,pos_y1:pos_y2) = ur_2;
        DEV.UR2(pos_x1:pos_x2,pos_y1:pos_y2) = ur; 
        DEV.UR2(:, pos_y2+1:Ny2) = ur2;

        DEV.RES = [dx dy];

        numgraph = 10;
        
        subplot(1,numgraph,1);

        xa = [-Nx2/2:Nx2/2]*dx2;
        ya = [0:Ny2-1]*dy2;
        [Y, X] = meshgrid(ya, xa);

        h = imagesc(xa,ya, DEV.UR2.',[1 10]);
        h2 = get(h, 'Parent');
        set(h2,'FontSize',10,'LineWidth',0.5);
        xlabel('$x$','Interpreter','LaTex');
        ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
        title('\mu_r');
        axis([-1 +1 -1 +1]);
        axis equal tight;
        colorbar;

        subplot(1,numgraph,2);

        g = imagesc(xa,ya, DEV.ER2.');
        g2 = get(g, 'Parent');
        set(g2,'FontSize',10,'LineWidth',0.5);
        xlabel('$x$','Interpreter','LaTex');
        ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
        title('\epsilon_r');
        axis([-1 +1 +0 +10])
        axis equal tight;
        colorbar;

        drawnow; 

        % The Two Dimensional Finite Difference Frequency Domain Method

        % NLAM = length(LAMBDA);

        THT = length(THETA);

        REF = zeros(THT,1);
        TRN = zeros(THT,1);
        CON = zeros(THT,1);

        for tht = 1:THT

            SRC.lam0 = LAMBDA(1);
            SRC.theta = THETA(tht);

            DAT = fdfd2d(DEV,SRC);
            % Result Translation

            Field = DAT.Field;
            R = DAT.RDE;
            T = DAT.TDE;
            REF(tht) = DAT.REF;
            TRN(tht) = DAT.TRN;
            CON(tht) = DAT.CON;
            

            % Displaying the Results

            %     index_m = Nx/2 + 1/2;
            %     disp('Reflection Diffraction Orders:' );
            %     disp(['RDE(-1) = '  num2str(100*R(index_m-1))]);
            %     disp(['RDE(0) = '  num2str(100*R(index_m))]);
            %     disp(['RDE(1) = '  num2str(100*R(index_m+1))]);
            %     disp(['RDE(2) = '  num2str(100*R(index_m+2))]);
            % 
            %     disp('Transmission Diffraction Orders:' );
            %     disp(['TDE(-1) = '  num2str(100*T(index_m-1))]);
            %     disp(['TDE(0) = '  num2str(100*T(index_m))]);
            %     disp(['TDE(1) = '  num2str(100*T(index_m+1))]);
            %     disp(['TDE(2) = '  num2str(100*T(index_m+2))]);


            disp(['REF(' num2str(tht) ') = ' num2str(REF(tht))]);
            disp(['TRN(' num2str(tht) ') = ' num2str(TRN(tht))]);
            disp(['CON(' num2str(tht) ') = ' num2str(CON(tht))]);

            % Post Processing Graphics

            subplot(1,numgraph,3);
            xb = [-floor(Nx/2):floor(Nx/2)]*dx;
            yb = [1:Ny]*dy;
            [Yb, Xb] = meshgrid(yb, xb);
            h = imagesc(xb,yb, real(Field).');
            h2 = get(h, 'Parent');
            set(h2,'FontSize',10,'LineWidth',0.5);
            xlabel('$x$','Interpreter','LaTex');
            ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
            mode = convertCharsToStrings(SRC.MODE);
            title([mode+ '-mode at ' + num2str(LAMBDA(1)/nanometers) + ' GHz', 'Re\{F\}']);
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
            title([mode+ '-mode at ' + num2str(LAMBDA(1)/nanometers) + ' GHz' , 'Im\{F\}']);
            shading interp;
            axis([-1 +1 -1 +1]);
            axis equal tight;
            colorbar;

            if THT <= 1 
                break; 
            end
            
            drawnow; %update graphics now!
           
        end
        
        REF2(index) = REF(tht);
        TRN2(index) = TRN(tht);
        CON2(index) = CON(tht);
        
    % For Convergence Test
        
        subplot(1,numgraph,5:10);
        plot(NRES(1:index),REF2(1:index),'r-o'); hold on;
        plot(NRES(1:index),TRN2(1:index),'b-o'); 
        plot(NRES(1:index),CON2(1:index),'k--o'); hold off;
        axis([min(NRES) max(NRES) 0 105]);
        xlabel('Grid Resolution');
        ylabel('Response ','Rotation',90);
        title('CONVERGENCE');
        
        h = legend('Reflectance','Transmittance', 'Conservation','Location','NorthEastOutside');
        set(h,'LineWidth',2);
        
        drawnow; %update graphics now!
        
    end

elseif FREQ == 1
    %% Frequency Code Block
    n_device = sqrt(ur*er); % refractive index of the device
    n_ref = sqrt(ur1*er1); % refractive index of the reflection region
    n_trn = sqrt(ur2*er2); % refractive index of the transmission region

    % Consider the Wavelength for resolving the resolution of the simulation

    n_max = max([n_device n_ref n_trn]);
    lambda_min = min(LAMBDA)/n_max;
    delta_gammaX = lambda_min/NRES(1);
    delta_gammaY = lambda_min/NRES(1);

    % Consider the Mechanical Features for optimizing the simulaiton
    mech_featX = [width L]; %array of mechanical features
    mech_featY = [height];
    d_minX = min(mech_featX)/4;
    d_minY = min(mech_featY)/4;

    % Choosing which is smaller between the parameters (to resolve smallest
    % features

    delta_x = min(delta_gammaX, d_minX);
    delta_y = min(delta_gammaY, d_minY);

    % Resolving Critical Dimensions through "Snapping"
    xTdev = min(mech_featX); %the critical dimension here is the duty cycle of the grating
    yTdev = min(mech_featY); %the height and the substrate thickness

    Mx = ceil(xTdev/delta_x);
    My = ceil(yTdev/delta_y);

    % Total Grid Size

    dx = xTdev/Mx; %the adjusted dx
    dy = yTdev/My; %the adjusted dy

    Nx = round(max(mech_featX)/dx)+ DEV.NPML(1)+ DEV.NPML(2)+ 2.*ceil((2*BUFZ)/(n_ref*dx));
    Nx = 2*round(Nx/2)+1; %physical size along x-axis
    Ny = round(sum(mech_featY)/dy)+ DEV.NPML(3)+ DEV.NPML(4)+ 2.*ceil((7*BUFZ)/(n_ref*dy)); %physical size along y-axis
    Ny = 2*round(Ny/2)+1;
    
    %The 2x Grid Parameters

    Nx2 = 2*Nx; %The 2x Grid requires 2x of the number of elements for both x and y
    Ny2 = 2*Ny;

    dx2 = dx/2;
    dy2 = dy/2;

    % Building the Simulation Matrix and the device
    
    pos_x0 = 2*DEV.NPML(1) + 2*round(BUFZ/dx2);
    pos_x0 = 2*round(pos_x0/2)+1;
    pos_x1 = pos_x0 + 1;
    pos_x2 = pos_x1 + round((L/dx2 - width/dx2)/2)-1;
    pos_x3 = pos_x2 + 1;
    pos_x4 = pos_x3 + round(width/dx2) - 1;
    pos_x5 = pos_x4 + 1;
    pos_x6 = pos_x5 + round((L/dx2 - width/dx2)/2)-1;
    pos_x7 = pos_x6 + 1;
    
    pos_y1 = DEV.NPML(3) + 3*round(BUFZ/dy2);
    pos_y1 = 2*round(pos_y1/2)+1;
    pos_y2 = pos_y1 + round(height/dy2) - 1;
    pos_y3 = pos_y2 + 1;
    pos_y4 = pos_y3 + round(subsheight/dy2) - 1; 
    
    DEV.ER2 = er1*ones(Nx2,Ny2);
    for index = 1:10:(pos_x7 - pos_x0) 
        %DEV.ER2(pos_x1+index:pos_x1+index+2,abs(round(pos_y1*(1-0.0000001*index.^2))):pos_y2) = er_2+0.015*index;
        DEV.ER2(pos_x1+index:pos_x1+index+2,abs(round(pos_y1*(1-0.00000001*index.^2))):pos_y2) = er_2;
        %DEV.ER2(pos_x1+index:pos_x1+index+2,pos_y3:abs(round(pos_y4*(1+0.0000001*index.^2)))) = er_2+0.010*index;
        DEV.ER2(pos_x1+index:pos_x7,pos_y3:pos_y4) = er_3;
    end
    %DEV.ER2(:, pos_y4+1:Ny2) = er2;
    
    DEV.UR2 = ur1*ones(Nx2,Ny2); 
%     for index = 1:(pos_x7 - pos_x0) 
%         DEV.UR2(pos_x1+index:pos_x7,pos_y1:pos_y2) = ur_2+0.010*index*rand(1,1);
%         DEV.UR2(pos_x1+index:pos_x7,pos_y3:pos_y4) = ur_2+0.015*index*rand(1,1);
%     end
%     DEV.UR2(:, pos_y4+1:Ny2) = ur2;
    
    DEV.RES = [dx dy];

    subplot(numgraph_row,numgraph_col,2);

    xa = [-Nx2/2:Nx2/2]*dx2;
    ya = [0:Ny2-1]*dy2;
    [Y, X] = meshgrid(ya, xa);

    h = imagesc(xa,ya, DEV.UR2.',[1 10]);
    h2 = get(h, 'Parent');
    set(h2,'FontSize',10,'LineWidth',0.5);
    xlabel('$x$','Interpreter','LaTex');
    ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
    title('\mu_r');
    axis([-1 +1 -1 +1]);
    axis equal tight;
    colorbar;

    subplot(numgraph_row,numgraph_col,1);

    g = imagesc(xa,ya, DEV.ER2.');
    g2 = get(g, 'Parent');
    set(g2,'FontSize',10,'LineWidth',0.5);
    xlabel('$x$','Interpreter','LaTex');
    ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
    title('\epsilon_r');
    axis([-1 +1 +0 +10])
    axis equal tight;
    colorbar;

    drawnow; 

% The Two Dimensional Finite Difference Frequency Domain Method

    NLAM = length(LAMBDA);
    REF = zeros(NLAM,1);
    TRN = zeros(NLAM,1);
    CON = zeros(NLAM,1);
    
    for nlam = 1:NLAM
        clearvars fd lamd x1 x2 x3 L d t ur er ur1 er1 ur2 er2 NRES BUFZ xbc ybc
        clearvars n_device n_ref n_trn n_max lambda_min
        clearvars delta_gammaX delta_gammaY mech_featX mech_featY
        clearvars d_minX d_minY delta_x delta_y xTdev yTdev
        clearvars Mx My dx2 dy2
        SRC.lam0 = LAMBDA(nlam);

        DAT = fdfd2d(DEV,SRC);
        
        clearvars ER2 UR2 NPML RES BC lam0 theta
        clearvars MODE Nx2 Ny2 k0 er_ref ur_ref er_trn ur_trn
        clearvars nref ntrn NGRID2X URxx URyy URzz 
        clearvars ERxx ERyy ERzz kinc m kx_m ky_ref_m
        clearvars ky_trn_m NGRID DEX DEY DHX DHY
        clearvars A xla yla Yla Xla f_src fsrc Q sr b Field
        clearvars Fref Ftrn x_value phaseTilt Aref Atrn
        clearvars Sref Strn
        
        % Result Translation

        Field = DAT.Field;
        R = DAT.RDE;
        T = DAT.TDE;
        REF(nlam) = DAT.REF;
        TRN(nlam) = DAT.TRN;
        CON(nlam) = DAT.CON;

        disp(['REF(' num2str(LAMBDA(nlam)/nanometers) ') = ' num2str(REF(nlam))]);
        disp(['TRN(' num2str(LAMBDA(nlam)/nanometers) ') = ' num2str(TRN(nlam))]);
        disp(['CON(' num2str(LAMBDA(nlam)/nanometers) ') = ' num2str(CON(nlam))]);

        % Post Processing Graphics

        subplot(numgraph_row,numgraph_col,3);
        xb = [-floor(Nx/2):floor(Nx/2)]*dx;
        yb = [1:Ny]*dy;
        [Yb, Xb] = meshgrid(yb, xb);
        h = imagesc(xb,yb, real(Field).');
        h2 = get(h, 'Parent');
        set(h2,'FontSize',10,'LineWidth',0.5);
        xlabel('$x$','Interpreter','LaTex');
        ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
        mode = convertCharsToStrings(SRC.MODE);
        title([mode+ '-mode at ' + num2str(LAMBDA(nlam)/nanometers) + 'nm', 'Re\{F\}']);
        axis([-1 +1 -1 +1]);
        axis equal tight;
        colorbar;
        %colormap(jet(1024));

        subplot(numgraph_row,numgraph_col,4);

        h = imagesc(xb,yb, imag(Field).');
        h2 = get(h, 'Parent');
        set(h2,'FontSize',10,'LineWidth',0.5);
        xlabel('$x$','Interpreter','LaTex');
        ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
        title([mode+ '-mode at ' + num2str(LAMBDA(nlam)/nanometers) + 'nm' , 'Im\{F\}']);
        shading interp;
        axis([-1 +1 -1 +1]);
        axis equal tight;
        colorbar;

        if NLAM <= 1
            break; 
        end
        if responsegraph == 1
            subplot(1,numgraph,5:10);

            plot(LAMBDA(1:nlam)./nanometers, REF(1:nlam),'-r'); hold on;
            plot(LAMBDA(1:nlam)./nanometers, TRN(1:nlam),'--b'); 
            plot(LAMBDA(1:nlam)./nanometers, CON(1:nlam),'-.k'); hold off;
            axis([min(LAMBDA/nanometers) max(LAMBDA/nanometers) 0 105]);
            xlabel('Wavelength (nm)');
            ylabel('Reflectance and Transmittance','Rotation',90);
            title('SPECTRAL RESPONSE')

            h = legend('Reflectance','Transmittance', 'Conservation','Location','NorthEastOutside');
            set(h,'LineWidth',2);
        end
        
        drawnow; %update graphics now!

    end
    
elseif ANGLE == 1
    %% Angle Sweep Code Block
    n_device = sqrt(ur*er); % refractive index of the device
    n_ref = sqrt(ur1*er1); % refractive index of the reflection region
    n_trn = sqrt(ur2*er2); % refractive index of the transmission region

    % Consider the Wavelength for resolving the resolution of the simulation

    n_max = max([n_device n_ref n_trn]);
    lambda_min = min(LAMBDA)/n_max;
    delta_gammaX = lambda_min/NRES(1);
    delta_gammaY = lambda_min/NRES(1);

    % Consider the Mechanical Features for optimizing the simulaiton
    mech_featX = [width L]; %array of mechanical features
    mech_featY = [height];
    d_minX = min(mech_featX)/4;
    d_minY = min(mech_featY)/4;

    % Choosing which is smaller between the parameters (to resolve smallest
    % features

    delta_x = min(delta_gammaX, d_minX);
    delta_y = min(delta_gammaY, d_minY);

    % Resolving Critical Dimensions through "Snapping"
    xTdev = min(mech_featX); %the critical dimension here is the duty cycle of the grating
    yTdev = min(mech_featY); %the height and the substrate thickness

    Mx = ceil(xTdev/delta_x);
    My = ceil(yTdev/delta_y);

    % Total Grid Size

    dx = xTdev/Mx; %the adjusted dx
    dy = yTdev/My; %the adjusted dy

    Nx = round(L/dx);
    Nx = 2*round(Nx/2)+1; %physical size along x-axis
    Ny = round(yTdev/dy)+ 2.*DEV.NPML(3)+ 2.*ceil(max(2*BUFZ)/(n_ref*dy)); %physical size along y-axis

    %The 2x Grid Parameters

    Nx2 = 2*Nx; %The 2x Grid requires 2x of the number of elements for both x and y
    Ny2 = 2*Ny;

    dx2 = dx/2;
    dy2 = dy/2;

    % Building the Simulation Matrix and the device

    pos_x1 = round((Nx2 - width/dx2)/2);
    pos_x2 = round((Nx2 + width/dx2)/2);
    
    pos_y1 = 2*DEV.NPML(3) + 2*round(BUFZ/dy2);
    pos_y1 = 2*round(pos_y1/2)+1;
    pos_y2 = pos_y1 + round(height/dy2) - 1;
    
    DEV.ER2 = er1*ones(Nx2,Ny2);
    DEV.ER2(:,pos_y1:pos_y2) = er_2;
    DEV.ER2(pos_x1:pos_x2,pos_y1:pos_y2) = er; 
    DEV.ER2(:, pos_y2+1:Ny2) = er2;
    
    DEV.UR2 = ur1*ones(Nx2,Ny2); %Since the relative permeability remains 1 (non-magnetic), fill the entire grid with ur1 = 1
    DEV.UR2(:,pos_y1:pos_y2) = ur_2;
    DEV.UR2(pos_x1:pos_x2,pos_y1:pos_y2) = ur; 
    DEV.UR2(:, pos_y2+1:Ny2) = ur2;
    
    DEV.RES = [dx dy];


    numgraph = 10;

    subplot(1,numgraph,1);

    xa = [-Nx2/2:Nx2/2]*dx2;
    ya = [0:Ny2-1]*dy2;
    [Y, X] = meshgrid(ya, xa);

    h = imagesc(xa,ya, DEV.UR2.',[1 10]);
    h2 = get(h, 'Parent');
    set(h2,'FontSize',10,'LineWidth',0.5);
    xlabel('$x$','Interpreter','LaTex');
    ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
    title('\mu_r');
    axis([-1 +1 -1 +1]);
    axis equal tight;
    colorbar;

    subplot(1,numgraph,2);

    g = imagesc(xa,ya, DEV.ER2.');
    g2 = get(g, 'Parent');
    set(g2,'FontSize',10,'LineWidth',0.5);
    xlabel('$x$','Interpreter','LaTex');
    ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
    title('\epsilon_r');
    axis([-1 +1 +0 +10])
    axis equal tight;
    colorbar;

    drawnow; 

    % The Two Dimensional Finite Difference Frequency Domain Method

    % NLAM = length(LAMBDA);

    THT = length(THETA);

    REF = zeros(THT,1);
    TRN = zeros(THT,1);
    CON = zeros(THT,1);

    for tht = 1:THT
        %clearvars fd lamd x1 x2 x3 L d t ur er ur1 er1 ur2 er2 NRES BUFZ xbc ybc
        clearvars n_device n_ref n_trn n_max lambda_min
        clearvars delta_gammaX delta_gammaY mech_featX mech_featY
        clearvars d_minX d_minY delta_x delta_y xTdev yTdev
        clearvars Mx My dx2 dy2
        SRC.lam0 = LAMBDA(1);
        SRC.theta = THETA(tht);

        DAT = fdfd2d(DEV,SRC);
        
        clearvars ER2 UR2 NPML RES BC lam0 theta
        clearvars MODE Nx2 Ny2 k0 er_ref ur_ref er_trn ur_trn
        clearvars nref ntrn NGRID2X sx sy URxx URyy URzz 
        clearvars ERxx ERyy ERzz kinc m kx_m ky_ref_m
        clearvars ky_trn_m NGRID DEX DEY DHX DHY
        clearvars A xla yla Yla Xla f_src fsrc Q sr b Field
        clearvars Fref Ftrn x_value phaseTilt Aref Atrn
        clearvars Sref Strn
        % Result Translation
        Field = DAT.Field;
        R = DAT.RDE;
        T = DAT.TDE;
        REF(tht) = DAT.REF;
        TRN(tht) = DAT.TRN;
        CON(tht) = DAT.CON;


        disp(['REF(' num2str(tht) ') = ' num2str(REF(tht))]);
        disp(['TRN(' num2str(tht) ') = ' num2str(TRN(tht))]);
        disp(['CON(' num2str(tht) ') = ' num2str(CON(tht))]);

        % Post Processing Graphics

        subplot(1,numgraph,3);
        xb = [-floor(Nx/2):floor(Nx/2)]*dx;
        yb = [1:Ny]*dy;
        [Yb, Xb] = meshgrid(yb, xb);
        h = imagesc(xb,yb, real(Field).');
        h2 = get(h, 'Parent');
        set(h2,'FontSize',10,'LineWidth',0.5);
        xlabel('$x$','Interpreter','LaTex');
        ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
        mode = convertCharsToStrings(SRC.MODE);
        title([mode+ '-mode at ' + num2str(LAMBDA(1)/nanometers) + ' nm', 'Re\{F\}']);
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
        title([mode+ '-mode at ' + num2str(LAMBDA(1)/nanometers) + ' nm' , 'Im\{F\}']);
        shading interp;
        axis([-1 +1 -1 +1]);
        axis equal tight;
        colorbar;

        if THT <= 1
            break; 
        end

        subplot(1,numgraph,5:10);
        
        plot(THETA(1:tht)/degrees, REF(1:tht),'-r'); hold on;
        plot(THETA(1:tht)/degrees, TRN(1:tht),'--b'); 
        plot(THETA(1:tht)/degrees, CON(1:tht),'-.k'); hold off;
        axis([min(THETA/degrees) max(THETA/degrees) 0 105]);
        xlabel('Angle of Incidence (deg)');
        ylabel('Reflectance and Transmittance','Rotation',90);
        title('SPECTRAL RESPONSE');
        
        h = legend('Reflectance','Transmittance', 'Conservation','Location','NorthEastOutside');
        set(h,'LineWidth',2);
        
        drawnow; %update graphics now!

    end
else
    disp('NO Simulation Procedure');
end

% Phase Animation and Movie Creation
% Draw Frame
if phasegraph == 1
    fig2 = figure('Color','w');
    set(fig2, 'Position', [1 41 1920 963]);
    set(fig2, 'Name', 'FDFD - Phase Animation');
    set(fig2, 'NumberTitle', 'off');
    set(fig2, 'Position', [1 41 1920 963]);
    % Draw Frames
    x_phasegraph = [-floor(Nx/2):floor(Nx/2)]*dx;
    y_phasegraph = [1:Ny]*dy;
    for nphi = 1:NPHI
        if MAKE_MOVIE == 1
            clf;
        end    
        field_phase = Field*exp(-1i*phase(nphi));
        imagesc(x_phasegraph, y_phasegraph, real(field_phase).');
        xlabel('$x$','Interpreter','LaTex');
        ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
        title('Steady-state Animation of the Field');
        shading interp;
        axis equal tight;
        colorbar;
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
end


% STOP TIMER
t2 = clock;
t = etime(t2,t1);
disp(['Elapsed time is ' num2str(t) ...
    ' seconds.']);


function DAT = fdfd2d(DEV,SRC)
% FDFD2D Two-Dimensional Finite-Difference Frequency-Domain
%
% DAT = fdfd2d(DEV,SRC)
%
% INPUT ARGUMENTS
% =================
% DEV Device Parameters
% .UR2 Relative permeability on 2X grid
% .ER2 Relative permittivity on 2X grid
% .NPML Size of PML on 1X grid [xlo xhi ylo yhi]
% .RES [dx2 dy2] grid resolution of 2X grid
%
% SRC
% .lam0 free space wavelength
% .theta Angle of incidence
% .MODE Mode: 'E' or 'H'
%
% OUTPUT ARGUMENTS
% =================
% DAT Output Data
% .RDE Array of diffraction efficiencies of reflected harmonics
% .REF Overall Reflectance
% .TDE Array of diffraction efficiencies of transmitted harmonics
% .TRN Overall Transmittance
% .CON Conservation of Energy
% .F Field
%
% Homework #8, Problem 1
% EE 5320 - COMPUTATIONAL ELECTROMAGNETICS

% Variable Translation
    ER2 = DEV.ER2;
    UR2 = DEV.UR2;
    NPML = DEV.NPML;
    RES = DEV.RES;
    BC = DEV.BC;
    scalar_pos = SRC.scalar_pos;

    lam0 = SRC.lam0;
    theta = SRC.theta;
    MODE = SRC.MODE;
    
    Nx = length(ER2(:,1))/2;
    Ny = length(ER2(1,:))/2;
    
    Nx2 = Nx*2;
    Ny2 = Ny*2;
    
    dx = RES(1);
    dy = RES(2);
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

    nref = sqrt(er_ref*ur_ref); %refractive index of the reflection region
    ntrn = sqrt(er_trn*ur_trn); %refractive index of the transmission region

    if er_ref <0 && ur_ref<0
        nref = - nref;
    end

    if er_trn <0 && ur_trn <0
        ntrn = - ntrn;
    end

% Calculate the Perfectly Matched Layer 

    NGRID2X = [Nx2 Ny2];
    [sx, sy] = calcpml2d(NGRID2X,NPML); % The function that calculates the Perfectly Matched Layer 
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
    if SRC.sourcetype == 0
        f_src = exp(1i.*(kinc(1).*X1a + kinc(2).*Y1a));
    end
    if SRC.sourcetype == 1
        [TH, R] = cart2pol(X1a,Y1a);
        [X1a, Y1a] = pol2cart(TH + theta, R);
        nxs = round(scalar_pos*Nx);
        X1a = X1a - X1a(nxs, 1);
        f_src = exp(-(2*X1a/SRC.w).^2).*exp(1i.*(kinc(2).*Y1a));%Source Field Matrix 
    end
    
    fsrc = f_src(:); %reshaping the Source Field Matrix into a column vector

% Compute the Scattered-Field Masking Matrix, Q

    Q = sparse(Nx, Ny); % Initiate the Q Matrix - Scattered-Field Masking Matrix
    sr = 1; % Set the amplitude to 1
    Q(:, 1:NPML(3)+2) = sr; %Fill the upper part of the 1x Grid with sr values
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
    DAT.Field = reshape(Field,Nx,Ny); %reshape the column vector to matrix form.

% POST PROCESSING: Computation of the Diffraction Efficiencies

% Extract Transmitted and Reflected Fields

    Fref = DAT.Field(:, NPML(3)+1); 
    Ftrn = DAT.Field(:, Ny - NPML(3)-1);

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
            DAT.RDE = abs(Sref).^2.*(real(ky_ref_m./ur_ref)./real(kinc(2)/ur_ref));
            DAT.TDE = abs(Strn).^2.*(real(ky_trn_m./ur_trn)./real(kinc(2)/ur_ref));
        case 'H'
            DAT.RDE = abs(Sref).^2.*(real(ky_ref_m./er_ref)./real(kinc(2)/er_ref));
            DAT.TDE = abs(Strn).^2.*(real(ky_trn_m./er_trn)./real(kinc(2)/er_ref));
    end

% Reflectance and Transmittance - Total of all the respective Diffraction
% Efficiencies

    DAT.REF = 100.*sum(DAT.RDE);
    DAT.TRN = 100.*sum(DAT.TDE);
    DAT.CON = DAT.REF + DAT.TRN;

end

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
        DEX = sparse(n,n);
        diagonal_0th = ones(n,1);
        diagonal_higher = ones(n,1);
        diagonal_higher(1:NGRID(1):NGRID(1)*NGRID(2)) = 0;
        DEX = spdiags([-diagonal_0th diagonal_higher], [0 1], DEX);
        
            
        if BC(1) == -2
            diagonal_lower = zeros(n,1);
            periodicityX = NGRID(1)*RES(1);
            diagonal_lower(1:NGRID(1):NGRID(1)*NGRID(2)) = exp(1i*kinc(1)*periodicityX);  
            DEX = spdiags([diagonal_lower], [-NGRID(1)+1], DEX);    
        end  
        DEX = DEX/RES(1);
    end

    if NGRID(2) == 1
        I = eye();
        DEY = 1i*kinc(2)*I;

    else
       n = NGRID(1)*NGRID(2);
       DEY = sparse(n,n);
       diagonal_0th = ones(n,1);
       diagonal_higher = ones(n,1);
       DEY = spdiags([-diagonal_0th diagonal_higher], [0 NGRID(1)], DEY);
            
       if BC(2) == -2
            diagonal_lower = zeros(n,1);
            periodicityY = NGRID(2)*RES(2);
            diagonal_lower(1:1:NGRID(1)*NGRID(2)) = exp(1i*kinc(2)*periodicityY);  
            DEY = spdiags([diagonal_lower], [-(NGRID(1)*NGRID(2)-(NGRID(1)-1))+1], DEY);
            
       end  
       DEY = DEY/RES(2);
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
        sx(NGRID(1)-NPML(2)+nx,:) = (1 + a_max*(nx/NPML(2))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*nx)/(2*NPML(2))))^2));

    end

    for ny = 1:NPML(3)
        sy(:,NPML(3)-ny+1) = (1 + a_max*(ny/NPML(3))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*ny)/(2*NPML(3))))^2));

    end

    for ny = 1 : NPML(4)
        sy(:,NGRID(2)-NPML(4)+ny) = (1 + a_max*(ny/NPML(4))^p)*(1 + 1i*eta0*(sigmaprime_max*(sin((pi*ny)/(2*NPML(4))))^2));

    end
end


