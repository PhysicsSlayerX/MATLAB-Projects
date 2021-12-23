% INITIALIZE MATLAB
close all;
clc;
clear all;

% TIME STAMP
t1 = clock;

% UNITS
degrees         = pi/180;
meters          = 1;
centimeters     = 1e-2 * meters;
millimmeters    = 1e-3 * meters;
micrometers     = 1e-6 * meters;
nanometers      = 1e-9 * meters;
inches          = 2.54 * centimeters;
feet            = 12 * inches;
seconds         = 1;
hertz           = 1/seconds;
kilohertz       = 1e3 * hertz;
megahertz       = 1e6 * hertz;
gigahertz       = 1e9 * hertz;
terahertz       = 1e12 * hertz;
petahertz       = 1e15 * hertz;

% CONSTANTS
e0 = 8.85418782e-12 * 1/meters;
u0 = 1.25663706e-6 * 1/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
NFREQ   = 100;
freq1   = 1* gigahertz;
freq2   = 5 * gigahertz;
FREQ    = linspace(freq1, freq2, NFREQ);

% LOSSY SLAB PARAMETERS

f0      = mean(FREQ);
lam0    = c0/f0;
% DIFFRACTION GRATING PARAMETERS
L       = 1.5 * centimeters;
%d       = 1.0 * centimeters;
d       = L;
er      = 15;
sig     = 150;

% PML PARAMETERS
pml_ky   = 1;
pml_ay   = 1e-10;
pml_Npml = 3;
pml_R0   = 1e-8;

% GRID PARAMETERS
lam_max = c0/min(FREQ);
NRES    = 132;
SPACER  = 0.05*lam_max * [1 1];
NPML    = [20 20];
nmax    = sqrt(er);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET FIGURE WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Color','k', 'Position', [0 42 1922 954]);
set(fig, 'Name', 'FDTD-2D Frequency Response Graph');
set(fig, 'NumberTitle', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTROLS FOR THE SHAPE GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_samples     = 1;

set_true            = 1;

plotControl         = 0;
imgImport           = 0;
centerRectangle     = set_true; %Shape1
centerCircle        = set_true; %Shape2
centerTriangle      = set_true; %Shape3
centerEllipse       = set_true; %Shape4
centerCross         = set_true; %Shape5
deformedCross       = set_true; %Shape6
rotatedTriangle     = set_true; %Shape7
rotatedRectangle    = set_true; %Shape8
rotatedEllipse      = set_true; %Shape9
combinedCircles     = set_true; %Shape10
rings               = set_true; %Shape11
centerTrap          = set_true; %Shape12
centerPent          = set_true; %Shape13
centerHex           = set_true; %Shape14
parallelogram_ring  = set_true; %Shape15
triangle_ring       = set_true; %Shape16
parallelogram_hole  = set_true; %Shape17
triangle_hole       = set_true; %Shape18
reg_hex_hole        = set_true; %Shape19
moon_shaped         = set_true; %Shape20
hex_circ_array      = set_true; %Shape21
binary_array        = 0; %Shape22
pixel_array         = 0; %Shape23

% Number of Samples per classification

class_sample    = 2500; % approx

step_var_rec    = ceil(class_sample^(1/2));       % Squared
step_var_circ   = ceil(class_sample);             % Single
step_var_tri    = ceil(class_sample^(1/3));       % Cubed
step_var_ell    = ceil(class_sample^(1/2));       % Squared
step_var_crs    = ceil(class_sample^(1/3));       % Cubed
step_var_dcrs   = ceil(class_sample^(1/3));       % Cubed
step_var_rtri   = ceil(class_sample^(1/2));       % Squared
step_var_rrect  = ceil(class_sample^(1/3));       % Cubed
step_var_rell   = ceil(class_sample^(1/3));       % Cubed
step_var_comc   = ceil(class_sample^(1/3));       % Cubed
step_var_ring   = ceil(class_sample^(1/3));       % Cubed
step_var_trap   = ceil(class_sample^(1/3));       % Cubed 
step_var_pent   = ceil(class_sample^(1/2));       % Squared
step_var_hex    = ceil(class_sample^(1/2));       % Squared
step_par_ring   = ceil(class_sample^(1/3));       % Cubed
step_tri_ring   = ceil(class_sample^(1/3));       % Cubed
step_par_hole   = ceil(class_sample^(1/3));       % Cubed
step_tri_hole   = ceil(class_sample^(1/3));       % Cubed
step_hex_hole   = ceil(class_sample^(1/3));       % Cubed
step_moon       = ceil(class_sample^(1/3));       % Cubed
step_hex_circ   = ceil(class_sample^(1/2));       % Squared
step_bin_array  = ceil(class_sample^(1/1));       % Single
step_pixel_arr  = ceil(class_sample^(1/1));       % Single

% Total Number of Samples
total_num_of_samples = centerRectangle*step_var_rec^2 + ...
                       centerCircle*step_var_circ + ...
                       centerTriangle*step_var_tri^3 + ...
                       centerEllipse*step_var_ell^2 +...
                       centerCross*step_var_crs^3 + ...
                       deformedCross*step_var_dcrs^3 + ...
                       rotatedTriangle*step_var_rtri^2 + ...
                       rotatedRectangle*step_var_rrect^3 + ...
                       rotatedEllipse*step_var_rell^3 + ...
                       combinedCircles*step_var_comc^3 + ...
                       rings*step_var_ring^3 + ...
                       centerTrap*step_var_trap^3 + ...
                       centerPent*step_var_pent^2 + ...
                       centerHex*step_var_hex^2 + ...
                       parallelogram_ring*step_par_ring^3 + ...
                       triangle_ring*step_tri_ring^3 + ...
                       parallelogram_hole*step_par_hole^3 + ...
                       triangle_hole*step_tri_hole^3 + ...
                       reg_hex_hole*step_hex_hole^3 + ...
                       moon_shaped*step_moon^3 + ...
                       hex_circ_array*step_hex_circ^2 + ...
                       binary_array*step_bin_array +...
                       pixel_array*step_pixel_arr;
                   

% Classification
num = 23;
classNum = linspace(1,num,num)-1;

% Seed Control
seedControl_trap     = 24;
seedControl_pent     = 2;
seedControl_hex      = 5;
seedControl_par_ring = 17;
seedControl_tri_ring = 96;
seedControl_par_hole = 4501;
seedControl_tri_hole = 450134;

% Display the Shape Collection Samples in a Group
example_width = 7;

% Number of elements for the binary elements (1,2,4,8,16,32,64,128)
num_elem_binary = 8;
num_elem_cont = 8;

% Put all the samples into one structure array
collate         = 1;

% Save the shapes
save_the_shapes = 0;
Filename = 'manyShapes_2500_per_type';
pixel_binary_group = 0;
Filename_pixel_binary = 'pixel_binary_figure';

% Save the shapes in HDF5 file format
save_the_shapes_HDF5 = 1;
filename_hdf5 = 'Shapes_2500_per_type_hdf5.h5';

% Save individual types

save_rec        = 0;       
save_circ       = 0;           
save_tri        = 0;      
save_ell        = 0;       
save_crs        = 0;       
save_dcrs       = 0;       
save_rtri       = 0;
save_rrect      = 0;
save_rell       = 0;
save_comc       = 0;
save_ring       = 0;
save_trap       = 0;
save_pent       = 0;
save_hex        = 0;
save_par_ring   = 0;
save_tri_ring   = 0;
save_par_hole   = 0;
save_tri_hole   = 0;
save_hex_hole   = 0;
save_moon       = 0;
save_hex_circ   = 0;
save_bin_array  = 0;
save_pixel_arr  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION
lam_min     = c0/max(FREQ);
dx          = lam_min/nmax/NRES;
dy          = lam_min/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);
dx = L/nx;
ny = ceil(d/dy);
dy = d/ny;

% COMPUTE GRID SIZE
Sx = L;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + d + SPACER(2);
Ny = NPML(1) + ceil(Sy/dy) + NPML(2);
Sy = Ny*dy;

% 2X GRID PARAMETERS
Nx2 = 2*Nx;     dx2 = dx/2;
Ny2 = 2*Ny;     dy2 = dy/2;

% GRID AXES
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ERzz    = ones(Nx,Ny);
URxx    = ones(Nx,Ny);
URyy    = ones(Nx,Ny);
SIGzz   = zeros(Nx,Ny);

f0  = 'DevER';              v0 = ones(Nx,Ny,1);
f01 = 'DevSIG';             v01 = zeros(Nx,Ny,1);
label = 'Label';            label01 = [];
image_num = 'Image_num';    image_num01 = [];

% Preallocate Spectral Structure Array (For Shapes)

Shape1 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape2 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape3 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape4 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape5 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape6 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape7 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape8 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape9 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape10 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape11 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape12 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape13 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape14 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape15 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape16 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape17 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape18 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape19 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape20 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape21 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape22 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);
Shape23 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMAGE IMPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if imgImport
    addpath 'D:\Research\Matlab practice\Numerical Transformation Optics\SNN_images'
    imagefiles = dir('D:\Research\Matlab practice\Numerical Transformation Optics\SNN_images\*.png');      
    nfiles = length(imagefiles);    % Number of files found
    for ii=1:nfiles
       currentfilename = imagefiles(ii).name;
       currentimage = imread(currentfilename);
       images{ii} = currentimage;
    end

    %nfiles = 5;
    for imageIndex = 1:nfiles
        %imageIndex = 400;
        thisImage = imrotate(images{imageIndex}, 90, 'bilinear', 'crop');
        thisImage = imresize(thisImage,[Nx Nx]);

        lastimage = ceil(thisImage(:,:,1));
        lastimage = im2double(lastimage);

        % BUILD IMAGE ARRAY
        ny1 = NPML(1) + round(0.75*SPACER(1)/dy) + 1;
        ny2 = ny1 + round(L/dx) - 1;

        % ADD PERMITTIVITY

        lastimage(lastimage > 0) = ceil(lastimage(lastimage > 0));
        ERzz(:,ny1:ny2) = er*lastimage;
        SIGzz(:,ny1:ny2) = sig*lastimage;

        ERzz = flip(ERzz, 1);
        SIGzz = flip(SIGzz, 1);
        ERzz(ERzz <= 1) = 1;
        SIGzz(SIGzz <= 0) = 0;

        %ERzz(ERzz > 1) = er;

        % SHOW DEVICE ERzz
        subplot(1, 2, 1);
        imagesc(xa, ya, ERzz.');
        axis equal tight off;
        title('ERzz');
        colorbar;
        plot_darkmode

        % SHOW DEVICE SIGzz
        subplot(1, 2, 2);
        imagesc(xa, ya, SIGzz.');
        axis equal tight off;
        title('SIGzz');
        colorbar;
        plot_darkmode;

        %drawnow;
        % EXTRACT ER and SIG MATRICES 
        imageER_zz      = mat2gray(ERzz(:,ny1:ny2));
        imageSIG_zz     = mat2gray(SIGzz(:,ny1:ny2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Centered Rectangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if centerRectangle
   

    % Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Customized number of samples (dual-looped)

    d1 = linspace(0.1,1.5,step_var_rec)*centimeters;

    for dimension2 = 1:length(d1)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);
        for dimension1 = 1:length(d1)

            % STOP/START INDICES
            nx = round(d1(dimension1)/dx);
            nx1 = 1 + floor((Nx - nx)/2);
            nx2 = nx1 + nx - 1;

            ny = round(d1(dimension2)/dy);
            ny1 = 1 + floor((Ny - ny)/2);
            ny2 = ny1 + ny - 1;


            % INCORPORATE RECTANGLE
            ERzz(nx1:nx2, ny1:ny2)  = er;
            SIGzz(nx1:nx2, ny1:ny2) = sig;

            if plotControl
                % SHOW DEVICE ERzz
                subplot(1, 2, 1);
                imagesc(xa, ya, ERzz.');
                axis equal tight off;
                title('ERzz');
                colorbar;
                plot_darkmode

                % SHOW DEVICE SIGzz
                subplot(1, 2, 2);
                imagesc(xa, ya, SIGzz.');
                axis equal tight off;
                title('SIGzz');
                colorbar;
                plot_darkmode;

                drawnow;
            end 

            imageER_zz      = ERzz(1:nx_design,ny_1d:ny_2d);
            imageSIG_zz     = SIGzz(1:nx_design,ny_1d:ny_2d);
            Shape1(dimension1+dimension2*step_var_rec - step_var_rec).DevER      = imageER_zz ;
            Shape1(dimension1+dimension2*step_var_rec - step_var_rec).DevSIG     = imageSIG_zz ;
            Shape1(dimension1+dimension2*step_var_rec - step_var_rec).Label      = classNum(1);
            Shape1(dimension1+dimension2*step_var_rec - step_var_rec).Image_num  = dimension1+dimension2*step_var_rec - step_var_rec;

            disp(['Processing Image: ' num2str(dimension1+dimension2*step_var_rec - step_var_rec) ' Class: ' num2str(classNum(1)) ])
        end

    end
    ShapeCollection = [Shape1];
    Centered_Rectangles = [Shape1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Centered Circles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if centerCircle
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % x and y axes

    xc = linspace(-L/2/centimeters, L/2/centimeters, Nx);
    yc = linspace(-L/2/centimeters, L/2/centimeters, Nx);

    % STOP/START INDICES
    nx = round(L/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;

    ny = round(L/dy);
    ny1 = 1 + floor((Ny - ny)/2);
    ny2 = ny1 + ny - 1;

    % INCORPORATE CIRCLE  

    [Xc, Yc] = meshgrid(xc, yc);

    Rc = linspace(0.1,L/2/centimeters,step_var_circ)*centimeters;

    for R_index = 1:length(Rc)
        ERzz(nx1:nx2, ny1:ny2) = (Xc.^2 + Yc.^2) <= (Rc(R_index)/centimeters)^2;
        ERzz(nx1:nx2, ny1:ny2) = er*ERzz(nx1:nx2, ny1:ny2);
        ERzz(ERzz <= 1) = 1;
        SIGzz(nx1:nx2, ny1:ny2) = (Xc.^2 + Yc.^2) <= (Rc(R_index)/centimeters)^2;
        SIGzz(nx1:nx2, ny1:ny2) = sig*SIGzz(nx1:nx2, ny1:ny2);
        SIGzz(SIGzz <= 0) = 0;

        if plotControl
            % SHOW DEVICE ERzz
            subplot(1, 2, 1);
            imagesc(xa, ya, ERzz.');
            axis equal tight off;
            title('ERzz');
            colorbar;
            plot_darkmode

            % SHOW DEVICE SIGzz
            subplot(1, 2, 2);
            imagesc(xa, ya, SIGzz.');
            axis equal tight off;
            title('SIGzz');
            colorbar;
            plot_darkmode;

            drawnow;
        end

        imageER_zz                 = ERzz(1:nx_design,ny_1d:ny_2d);
        imageSIG_zz                = SIGzz(1:nx_design,ny_1d:ny_2d);
        Shape2(R_index).DevER      = imageER_zz ;
        Shape2(R_index).DevSIG     = imageSIG_zz ;
        Shape2(R_index).Label      = classNum(2);
        Shape2(R_index).Image_num  = step_var_rec^2 + R_index;

        disp(['Processing Image: ' num2str(R_index) ' Radius: ' num2str(Rc(R_index)/centimeters) ' cm ' 'Class: ' num2str(classNum(2))])
    end
    ShapeCollection = [Shape2];
    Centered_Circles = [Shape2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Centered Triangles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if centerTriangle

    % Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;

    tri_1 = linspace(0.3,1.5,step_var_tri)*centimeters;

    % dimension factors
    width_factor    = linspace(0.3,1.0,step_var_tri);
    height_factor   = linspace(0.3,1.0,step_var_tri);

    % INCORPORATE TRIANGLE
    for h_index = 1:length(height_factor)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);
        for w_index = 1:length(width_factor)
            ERzz    = ones(Nx,Ny);
            URxx    = ones(Nx,Ny);
            URyy    = ones(Nx,Ny);
            SIGzz   = zeros(Nx,Ny);
            for t_index = 1:length(tri_1)

                % STOP/START INDICES
                h       = tri_1(t_index)*sqrt(3)/2;
                ny      = height_factor(h_index)*round(h/dy);
                ny1     = 1 + floor((Ny - ny)/2);
                ny2     = ny1 + ny - 1;

                for ny = ny1 : ny2

                    % Calculate the 0 to 1 number
                    f = width_factor(w_index)*(ny - ny1 + 1)/(ny2 - ny1 + 1);

                    % Calculate Width of Triangle at Current Position
                    nx   = round(f*(tri_1(t_index))/dx);
                    nx_d = round(f*(tri_1(t_index))/dx);

                    % Calculate Stop/Start Indices
                    nx1 = 1 + floor((Nx - nx)/2);
                    nx2 = nx1 + nx - 1;

                    % Incorporate 1's into A
                    ERzz(nx1:nx2,ny)    = er;
                    SIGzz(nx1:nx2,ny)   = sig;
                end

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                    
                end 
                
                indexName                    = (t_index + w_index*step_var_tri - step_var_tri) + h_index*step_var_tri^2 - step_var_tri^2;
                
                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape3(indexName).DevER      = imageER_zz ;
                Shape3(indexName).DevSIG     = imageSIG_zz ;
                Shape3(indexName).Label      = classNum(3);
                
                indexLabel                   = step_var_rec^2 + step_var_circ + indexName;
                
                Shape3(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(3)) ' Label ' num2str(indexLabel)])
           
            end
        end
    end
    ShapeCollection = [Shape3];
    Centered_Triangles = [Shape3];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Centered Ellipses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if centerEllipse
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;

    Shape4 = struct(f0,v0,f01,v01, label, label01, image_num, image_num01);

    xc = linspace(-L/2/centimeters, L/2/centimeters, Nx);
    yc = linspace(-L/2/centimeters, L/2/centimeters, Nx);

    % STOP/START INDICES
    nx = round(L/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;

    ny = round(L/dy);
    ny1 = 1 + floor((Ny - ny)/2);
    ny2 = ny1 + ny - 1;

    % INCORPORATE CIRCLE  

    [Xc, Yc] = meshgrid(xc, yc);
    
    rx = linspace(0.1,L/2/centimeters,step_var_ell);
    ry = linspace(0.1,L/2/centimeters,step_var_ell);
    Rc = L/2/centimeters; %Fixed this value to be the radius of the largest circle that will fit inside the square.

    for ry_index = 1:length(ry)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);

        for rx_index = 1:length(rx)
            
            if ry_index ~= rx_index
                ERzz(nx1:nx2, ny1:ny2) = ((Xc/rx(rx_index)).^2 + (Yc/ry(ry_index)).^2) <= (Rc)^2;
                ERzz(nx1:nx2, ny1:ny2) = er*ERzz(nx1:nx2, ny1:ny2);
                ERzz(ERzz <= 1) = 1;
                SIGzz(nx1:nx2, ny1:ny2) = ((Xc/rx(rx_index)).^2 + (Yc/ry(ry_index)).^2) <= (Rc)^2;
                SIGzz(nx1:nx2, ny1:ny2) = sig*SIGzz(nx1:nx2, ny1:ny2);
                SIGzz(SIGzz <= 0) = 0;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end

                
            end
            indexName                  = rx_index + ry_index*step_var_ell - step_var_ell;

            imageER_zz                 = ERzz(1:nx_design,ny_1d:ny_2d);
            imageSIG_zz                = SIGzz(1:nx_design,ny_1d:ny_2d);
            Shape4(indexName).DevER      = imageER_zz ;
            Shape4(indexName).DevSIG     = imageSIG_zz ;
            Shape4(indexName).Label      = classNum(4);

            indexLabel                 = step_var_circ + step_var_rec^2 + step_var_tri^3 + indexName;
            Shape4(indexName).Image_num  = indexLabel; 

            disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(4)) ' Label ' num2str(indexLabel)])
        end
    end
    ShapeCollection = [Shape4];
    Centered_Ellipses = [Shape4];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Centered Crosses(Uniform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if centerCross
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    

    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Changing Values of nx_v, ny_v, and phi
    leg1 = linspace(L/15,L/2,step_var_crs);
    leg2 = linspace(L/15,L/2,step_var_crs);
    rot_crs = linspace(0,90,step_var_crs)*degrees;
    
    
    for ldim1 = 1:length(leg1)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);
        
        for ldim2 = 1:length(leg2)
            
            for rot_dim = 1:length(rot_crs)
            
                [TY, TX] = meshgrid(ya1, xa1);

                % FIRST PART
                % ROTATE COORDINATES

                phi1             = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES
                nx_v1 = leg1(ldim1);        % handles
                nx1_v1 = -nx_v1;
                nx2_v1 = nx_v1;

                ny_v1 = leg2(ldim2);        % handles
                ny1_v1 = -ny_v1;
                ny2_v1 = ny_v1;

                TempArray1 = (TX >= nx1_v1 & TX<= nx2_v1 & TY >= ny1_v1 & TY <= ny2_v1);

                % SECOND PART
                % ROTATE COORDINATES
                %[TY, TX] = meshgrid(ya1, xa1);

                phi2            = 90*degrees;   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi2, RADIUS);

                % PRECISE VERTICES
                nx_v2 = nx_v1;        % handles
                nx1_v2 = -nx_v2;
                nx2_v2 = nx_v2;

                ny_v2 = ny_v1;        % handles
                ny1_v2 = -ny_v2;
                ny2_v2 = ny_v2;

                TempArray2 = (TX >= nx1_v2 & TX<= nx2_v2 & TY >= ny1_v2 & TY <= ny2_v2);

                TempArray = TempArray1 | TempArray2;

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*double(TempArray);
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*double(TempArray);
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim2*step_var_crs - step_var_crs) + ldim1*step_var_crs^2 - step_var_crs^2;
                
                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape5(indexName).DevER      = imageER_zz ;
                Shape5(indexName).DevSIG     = imageSIG_zz ;
                Shape5(indexName).Label      = classNum(5);
                
                indexLabel                   = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 + indexName;
                
                Shape5(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(5)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape5];
    Centered_Crosses = [Shape5];
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Deformed Crosses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if deformedCross
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    
    % Changing Values of nx_v, ny_v, and phi
    leg1 = linspace(L/15,L/2,step_var_dcrs);
    leg2 = linspace(L/15,L/2,step_var_dcrs);
    rot_crs = linspace(0,180,step_var_dcrs)*degrees;

    
    for ldim1 = 1:length(leg1)        
        for ldim2 = 1:length(leg2)       
            for rot_dim = 1:length(rot_crs)
            
                [TY, TX] = meshgrid(ya1, xa1);

                % FIRST PART
                % ROTATE COORDINATES

                phi1             = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES
                nx_v1 = leg1(ldim1);        % handles
                nx1_v1 = -nx_v1;
                nx2_v1 = nx_v1;

                ny_v1 = leg2(ldim2);        % handles
                ny1_v1 = -ny_v1;
                ny2_v1 = ny_v1;

                TempArray1 = (TX >= nx1_v1 & TX<= nx2_v1 & TY >= ny1_v1 & TY <= ny2_v1);

                % SECOND PART
                % ROTATE COORDINATES
                %[TY, TX] = meshgrid(ya1, xa1);

                phi2            = rot_crs(rot_dim)+90*degrees;   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi2, RADIUS);

                % PRECISE VERTICES
                nx_v2 = nx_v1;        % handles
                nx1_v2 = -nx_v2;
                nx2_v2 = nx_v2;

                ny_v2 = ny_v1;        % handles
                ny1_v2 = -ny_v2;
                ny2_v2 = ny_v2;

                TempArray2 = (TX >= nx1_v2 & TX<= nx2_v2 & TY >= ny1_v2 & TY <= ny2_v2);

                TempArray = TempArray1 | TempArray2;

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*double(TempArray);
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*double(TempArray);
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim2*step_var_dcrs - step_var_dcrs) + ldim1*step_var_dcrs^2 - step_var_dcrs^2;
                
                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape6(indexName).DevER      = imageER_zz ;
                Shape6(indexName).DevSIG     = imageSIG_zz ;
                Shape6(indexName).Label      = classNum(6);
                
                indexLabel                   = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 + step_var_crs^3 + indexName;
                
                Shape6(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(6)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape6];
    Deformed_Crosses = [Shape6];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Rotated Triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rotatedTriangle
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);

    % Handle Values
    RAD = linspace(L/10,L/2,step_var_rtri);
    rot_crs = linspace(0,90,step_var_rtri)*degrees;
    rand_rot = randi([115, 125], 1, step_var_rtri);

    
    for ldim1 = 1:length(RAD)             
        for rot_dim = 1:length(rot_crs)
            [TY, TX] = meshgrid(ya1, xa1);
            % ROTATE COORDINATES

            phi1            = rot_crs(rot_dim);   % handles
            [THETA, RADIUS] = cart2pol(TX, TY);
            [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

            % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

            %RAD = L/3; %Distance from (0,0);

            p1 = [(RAD(ldim1))*cos(phi1) ; (RAD(ldim1))*sin(phi1)];
            p2 = [RAD(ldim1)*cos(phi1+rand_rot(rot_dim)*degrees) ; RAD(ldim1)*sin(phi1+rand_rot(rot_dim)*degrees)];
            p3 = [RAD(ldim1)*cos(phi1+2*rand_rot(rot_dim)*degrees) ; RAD(ldim1)*sin(phi1+2*rand_rot(rot_dim)*degrees)];

            P  = [ p1 p2 p3 ];

            A = polyfill(xa1,ya1,P);
            A = fliplr(A); % use fliplr rather tha flipud

            % STOP/START INDICES
            nx = round(L/dx);
            nx1 = 1 + floor((Nx - nx)/2);
            nx2 = nx1 + nx - 1;

            ny = round(L/dy);
            ny1 = 1 + floor((Ny - ny)/2);
            ny2 = ny1 + ny - 1;

            % INCORPORATE MATERIAL PROPERTIES
            ERzz(nx1:nx2, ny1:ny2)  = er*A;
            ERzz(ERzz <= 1)         = 1;
            SIGzz(nx1:nx2, ny1:ny2) = sig*A;
            SIGzz(SIGzz <= 1)       = 1;

            if plotControl
                % SHOW DEVICE ERzz
                subplot(1, 2, 1);
                imagesc(xa, ya, ERzz.');
                axis equal tight off;
                title('ERzz');
                colorbar;
                plot_darkmode

                % SHOW DEVICE SIGzz
                subplot(1, 2, 2);
                imagesc(xa, ya, SIGzz.');
                axis equal tight off;
                title('SIGzz');
                colorbar;
                plot_darkmode;

                drawnow;
            end 

            indexName                    = (rot_dim + ldim1*step_var_rtri - step_var_rtri);

            imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
            imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
            Shape7(indexName).DevER      = imageER_zz ;
            Shape7(indexName).DevSIG     = imageSIG_zz ;
            Shape7(indexName).Label      = classNum(7);

            indexLabel                   = step_var_rec^2 + step_var_circ + step_var_tri^3 ...
                                           + step_var_ell^2 + step_var_crs^3 + ...
                                           + step_var_dcrs^3 + indexName;

            Shape7(indexName).Image_num  = indexLabel;

            disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(7)) ' Label ' num2str(indexLabel)])
        end
    end
    ShapeCollection = [Shape7];
    Rotated_Triangles = [Shape7];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Rotated Rectangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rotatedRectangle
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);

    % Handle Values
    rrectdim1 = linspace(L/10,L/3.5,step_var_rrect);
    rrectdim2 = linspace(L/10,L/3.5,step_var_rrect);
    rot_crs = linspace(0,180,step_var_rrect)*degrees;
   

    for ldim1 = 1:length(rrectdim1)     
        for ldim2 = 1:length(rrectdim2)  
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);

                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES
                nx_v1 = rrectdim1(ldim1);        % handles
                nx1_v1 = -nx_v1;
                nx2_v1 = nx_v1;

                ny_v1 = rrectdim2(ldim2);        % handles
                ny1_v1 = -ny_v1;
                ny2_v1 = ny_v1;

                TmpArr = (TX >= nx1_v1 & TX<= nx2_v1 & TY >= ny1_v1 & TY <= ny2_v1);

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*TmpArr;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*TmpArr;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim2*step_var_rrect - step_var_rrect) + ldim1*step_var_rrect^2 - step_var_rrect^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape8(indexName).DevER      = imageER_zz ;
                Shape8(indexName).DevSIG     = imageSIG_zz ;
                Shape8(indexName).Label      = classNum(8);

                indexLabel                   = step_var_rec^2 + step_var_circ + step_var_tri^3 ...
                                               + step_var_ell^2 + step_var_crs^3 + ...
                                               + step_var_dcrs^3 + step_var_rtri^2 + indexName;

                Shape8(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(8)) ' Label ' num2str(indexLabel)])
                
            end
        end
    end
    ShapeCollection = [Shape8];
    Rotated_Rectangle = [Shape8];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Rotated Ellipses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rotatedEllipse
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    xc = [0:Nx-1]*dx;   xc = xc - mean(xc);
    yc = [0:Nx-1]*dy;   yc = yc - mean(yc);
    
    % STOP/START INDICES
    nx = round(L/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;

    ny = round(L/dy);
    ny1 = 1 + floor((Ny - ny)/2);
    ny2 = ny1 + ny - 1;

    rx = linspace(0.1,L/2/centimeters,step_var_rell);
    ry = linspace(0.1,L/2/centimeters,step_var_rell);
    Rc = L/2/centimeters; %Fixed this value to be the radius of the largest circle that will fit inside the square.
    rot_crs = linspace(0,180,step_var_rell)*degrees;
    
    for rx_index = 1:length(rx)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);

        for ry_index = 1:length(ry)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);
               
            for rot_dim = 1:length(rot_crs)

                % INCORPORATE ROTATING ELLIPSE  
                [Yc, Xc] = meshgrid(yc, xc);

                phi_ellp        = rot_crs(rot_dim);
                [THETA, RADIUS] = cart2pol(Xc, Yc);
                [Xc, Yc]        = pol2cart(THETA+phi_ellp, RADIUS);

                Xc = Xc./centimeters;
                Yc = Yc./centimeters;

                ERzz(nx1:nx2, ny1:ny2) = ((Xc/rx(rx_index)).^2 + (Yc/ry(ry_index)).^2) <= (Rc)^2;
                ERzz(nx1:nx2, ny1:ny2) = er*ERzz(nx1:nx2, ny1:ny2);
                ERzz(ERzz <= 1) = 1;
                SIGzz(nx1:nx2, ny1:ny2) = ((Xc/rx(rx_index)).^2 + (Yc/ry(ry_index)).^2) <= (Rc)^2;
                SIGzz(nx1:nx2, ny1:ny2) = sig*SIGzz(nx1:nx2, ny1:ny2);
                SIGzz(SIGzz <= 0) = 0;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end
               
                % Display in Command Window and Save the Image in the
                % corresponding structure array
                indexName                  = (rot_dim + ry_index*step_var_rell - step_var_rell) + rx_index*step_var_rell^2 - step_var_rell^2;

                imageER_zz                 = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape9(indexName).DevER      = imageER_zz ;
                Shape9(indexName).DevSIG     = imageSIG_zz ;
                Shape9(indexName).Label      = classNum(9);

                indexLabel                 = step_var_rec^2 + step_var_circ + step_var_tri^3 ...
                                             + step_var_ell^2 + step_var_crs^3 + ...
                                             + step_var_dcrs^3 + step_var_rtri^2 + ...
                                             + step_var_rrect^3 + indexName;
                Shape9(indexName).Image_num  = indexLabel; 

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(9)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape9];
    Rotated_Ellipses = [Shape9];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Combined Circles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if combinedCircles
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;

    xc = [0:Nx-1]*dx;   xc = xc - mean(xc);
    yc = [0:Nx-1]*dy;   yc = yc - mean(yc);
    
    % STOP/START INDICES
    nx = round(L/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;

    ny = round(L/dy);
    ny1 = 1 + floor((Ny - ny)/2);
    ny2 = ny1 + ny - 1;

    rx = linspace(-0.1,-L/3/centimeters,step_var_comc);
    Rc = linspace(0.1,L/3/centimeters,step_var_comc);
    rot_crs = linspace(0,180,step_var_comc)*degrees;
    
    %rot_crs = 0*degrees;
    
    for rx_index = 1:length(rx)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);

        for rc_index = 1:length(Rc)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);
               
            for rot_dim = 1:length(rot_crs)

                % INCORPORATE ROTATING ELLIPSE  
                [Yc, Xc] = meshgrid(yc, xc);

                phi_comc        = rot_crs(rot_dim);
                [THETA, RADIUS] = cart2pol(Xc, Yc);
                [Xc, Yc]        = pol2cart(THETA+phi_comc, RADIUS);

                Xc = Xc./centimeters;
                Yc = Yc./centimeters;
                
                circ1 = ((Xc+rx(rx_index)).^2 + (Yc).^2) <= (Rc(rc_index))^2;
                circ2 = ((Xc-rx(rx_index)).^2 + (Yc).^2) <= (Rc(rc_index))^2;
             

                ERzz(nx1:nx2, ny1:ny2) = double(circ1 | circ2);
                ERzz(nx1:nx2, ny1:ny2) = er*ERzz(nx1:nx2, ny1:ny2);
                ERzz(ERzz <= 1) = 1;
                SIGzz(nx1:nx2, ny1:ny2) = double(circ1 | circ2);
                SIGzz(nx1:nx2, ny1:ny2) = sig*SIGzz(nx1:nx2, ny1:ny2);
                SIGzz(SIGzz <= 0) = 0;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end
               
                % Display in Command Window and Save the Image in the
                % corresponding structure array
                indexName                  = (rot_dim + rc_index*step_var_comc - step_var_comc) + rx_index*step_var_comc^2 - step_var_comc^2;

                imageER_zz                 = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape10(indexName).DevER      = imageER_zz ;
                Shape10(indexName).DevSIG     = imageSIG_zz ;
                Shape10(indexName).Label      = classNum(10);

                indexLabel                 = step_var_rec^2 + step_var_circ + step_var_tri^3 ...
                                             + step_var_ell^2 + step_var_crs^3 + ...
                                             + step_var_dcrs^3 + step_var_rtri^2 + ...
                                             + step_var_rrect^3 + step_var_rell^3 + indexName;
                Shape10(indexName).Image_num  = indexLabel; 

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(10)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape10];
    Combined_Circles = [Shape10];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Elliptical and Circular Rings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rings
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    %Design Space
    nx_design = round(L/dx);

    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    xc = [0:Nx-1]*dx;   xc = xc - mean(xc);
    yc = [0:Nx-1]*dy;   yc = yc - mean(yc);
    
    %STOP/START INDICES
    nx = round(L/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;

    ny = round(L/dy);
    ny1 = 1 + floor((Ny - ny)/2);
    ny2 = ny1 + ny - 1;
    
    rx       = linspace(0.1,L/2/centimeters,step_var_ring);
    ry       = linspace(0.1,L/2/centimeters,step_var_ring); 
    rx2      = linspace(0.1,L/2/centimeters,step_var_ring);
    ry2      = linspace(0.1,L/2/centimeters,step_var_ring);
    Rc       = L/2/centimeters;
    scale_f  = 0.5;
    rot_crs  = linspace(0,180,step_var_ring)*degrees;
    
    for rx_index = 1:length(rx)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);

        for ry_index = 1:length(ry)
        ERzz    = ones(Nx,Ny);
        URxx    = ones(Nx,Ny);
        URyy    = ones(Nx,Ny);
        SIGzz   = zeros(Nx,Ny);
               
            for rot_dim = 1:length(rot_crs)

                %INCORPORATE ROTATING ELLIPSE  
                [Yc, Xc] = meshgrid(yc, xc);

                phi_comc        = rot_crs(rot_dim);
                [THETA, RADIUS] = cart2pol(Xc, Yc);
                [Xc, Yc]        = pol2cart(THETA+phi_comc, RADIUS);

                Xc = Xc./centimeters;
                Yc = Yc./centimeters;
                
                ellipse1 = (Xc/rx(rx_index)).^2 + (Yc/ry(ry_index)).^2 <= Rc.^2;
                ellipse2 = (Xc/(rx2(rx_index)*scale_f)).^2 + (Yc/(ry2(ry_index)*scale_f)).^2 <= Rc.^2;  
             

                ERzz(nx1:nx2, ny1:ny2) = double(xor(ellipse1,ellipse2));
                ERzz(nx1:nx2, ny1:ny2) = er*ERzz(nx1:nx2, ny1:ny2);
                ERzz(ERzz <= 1) = 1;
                SIGzz(nx1:nx2, ny1:ny2) = double(xor(ellipse1,ellipse2));
                SIGzz(nx1:nx2, ny1:ny2) = sig*SIGzz(nx1:nx2, ny1:ny2);
                SIGzz(SIGzz <= 0) = 0;

                if plotControl
                    %SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    %SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end
               
                %Display in Command Window and Save the Image in the
                %corresponding structure array
                indexName                     = (rot_dim + ry_index*step_var_ring - step_var_ring) + rx_index*step_var_ring^2 - step_var_ring^2;

                imageER_zz                    = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                   = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape11(indexName).DevER      = imageER_zz ;
                Shape11(indexName).DevSIG     = imageSIG_zz ;
                Shape11(indexName).Label      = classNum(11);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + indexName;
                Shape11(indexName).Image_num  = indexLabel; 

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(11)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape11];
    Elliptic_Rings  = [Shape11];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Centered Trapezoids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if centerTrap
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);

    % Handle Values
    
    limitA = L/15;
    limitB = L/3;
    N      = step_var_trap;
    rng(seedControl_trap, 'combRecursive');
    
    rand1x = limitA + (limitB - limitA).*rand(1,N);
    rand1y = limitA + (limitB - limitA).*rand(1,N);
    rand2x = limitA + (limitB - limitA).*rand(1,N);
    rand2y = limitA + (limitB - limitA).*rand(1,N);
    rand3x = limitA + (limitB - limitA).*rand(1,N);
    rand3y = limitA + (limitB - limitA).*rand(1,N);
    rand4x = limitA + (limitB - limitA).*rand(1,N);
    rand4y = limitA + (limitB - limitA).*rand(1,N);
    
    rot_crs = linspace(0,180,step_var_trap)*degrees;
    
    for ldim1 = 1:length(rand1x)     
        for ldim2 = 1:length(rand1x)  
            for rot_dim = 1:length(rot_crs)
                
                [TY, TX] = meshgrid(ya1, xa1);

                % ROTATE COORDINATES
                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % DEFINE TWO POINTS
                p1 = [-rand1x(ldim1), +rand1y(ldim2)];
                p2 = [+rand2x(ldim1), +rand2y(ldim2)];
                p3 = [-rand3x(ldim1), -rand3y(ldim2)];
                p4 = [+rand4x(ldim1), -rand4y(ldim2)];
               
                
                m1 = (p2(2)-p1(2))/(p2(1) - p1(1));
                m2 = (p2(2)-p4(2))/(p2(1) - p4(1));
                m3 = (p3(2)-p4(2))/(p3(1) - p4(1));
                m4 = (p3(2)-p1(2))/(p3(1) - p1(1));
                
                % FILL HALF SPACES
                TmpArr1 = (TY - p1(2)) - m1*(TX  - p1(1)) < 0 ;
                if p2(1) < p4(1)
                    TmpArr2 = (TY - p4(2)) - m2*(TX  - p4(1)) < 0;
                elseif p2(1) == p4(1)
                    TmpArr2 = (TY - p4(2)) - m2*(TX  - p4(1)) > 0;
                else
                    TmpArr2 = (TY - p4(2)) - m2*(TX  - p4(1)) > 0;
                end
                TmpArr3 = (TY - p4(2)) - m3*(TX  - p4(1)) > 0;
                if p1(1) <= p3(1)
                    TmpArr4 = (TY - p1(2)) - m4*(TX  - p1(1)) > 0;
                else
                    TmpArr4 = (TY - p1(2)) - m4*(TX  - p1(1)) < 0;
                end
                
                TmpArr = (TmpArr1 .* TmpArr2 .* TmpArr3 .* TmpArr4); 
                TmpArr = fliplr(TmpArr);
                
                trap_plot_test = 0;
                
                if trap_plot_test
                    % SHOW DEVICE ERzz
                    subplot(1, 5, 1);
                    imagesc(xa1, ya1, fliplr(TmpArr1).');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    subplot(1, 5, 2);
                    imagesc(xa1, ya1, fliplr(TmpArr2).');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    subplot(1, 5, 3);
                    imagesc(xa1, ya1, fliplr(TmpArr3).');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    subplot(1, 5, 4);
                    imagesc(xa1, ya1, fliplr(TmpArr4).');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    subplot(1, 5, 5);
                    imagesc(xa1, ya1, TmpArr.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode
                    
                    drawnow;
                end
                
                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*TmpArr;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*TmpArr;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim2*step_var_trap - step_var_trap) + ldim1*step_var_trap^2 - step_var_trap^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape12(indexName).DevER      = imageER_zz ;
                Shape12(indexName).DevSIG     = imageSIG_zz ;
                Shape12(indexName).Label      = classNum(12);

                indexLabel                   = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                               step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                               step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + indexName;

                Shape12(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(12)) ' Label ' num2str(indexLabel)])
                
            end
        end
    end
    rng default
    ShapeCollection = [Shape12];
    Centered_Trapezoids = [Shape12];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Centered Pentagon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if centerPent
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Random Values for R
    
    limitA = L/10;
    limitB = L/2;
    rng(seedControl_pent);
    RAD1x = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD2x = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD3x = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD4x = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD5x = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD1y = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD2y = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD3y = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD4y = limitA + (limitB - limitA).*rand(1,step_var_pent);
    RAD5y = limitA + (limitB - limitA).*rand(1,step_var_pent);

    
    % Handle Values
    rot_crs = linspace(0,90,step_var_pent)*degrees;

    for ldim1 = 1:length(RAD1x)             
        for rot_dim = 1:length(rot_crs)
            [TY, TX] = meshgrid(ya1, xa1);
            % ROTATE COORDINATES

            phi1            = rot_crs(rot_dim);   % handles
            [THETA, RADIUS] = cart2pol(TX, TY);
            [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

            % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

            %RAD = L/3; %Distance from (0,0);

            p1 = [RAD1x(ldim1)*cos(phi1)              ; RAD1y(ldim1)*sin(phi1)];
            p2 = [RAD2x(ldim1)*cos(phi1+72*degrees)   ; RAD2y(ldim1)*sin(phi1+72*degrees)];
            p3 = [RAD3x(ldim1)*cos(phi1+72*2*degrees) ; RAD3y(ldim1)*sin(phi1+72*2*degrees)];
            p4 = [RAD4x(ldim1)*cos(phi1+72*3*degrees) ; RAD4y(ldim1)*sin(phi1+72*3*degrees)];
            p5 = [RAD5x(ldim1)*cos(phi1+72*4*degrees) ; RAD5y(ldim1)*sin(phi1+72*4*degrees)];
            
            P  = [ p1 p2 p3 p4 p5];

            A = polyfill(xa1,ya1,P);
            A = fliplr(A); % use fliplr rather tha flipud

            % STOP/START INDICES
            nx = round(L/dx);
            nx1 = 1 + floor((Nx - nx)/2);
            nx2 = nx1 + nx - 1;

            ny = round(L/dy);
            ny1 = 1 + floor((Ny - ny)/2);
            ny2 = ny1 + ny - 1;

            % INCORPORATE MATERIAL PROPERTIES
            ERzz(nx1:nx2, ny1:ny2)  = er*A;
            ERzz(ERzz <= 1)         = 1;
            SIGzz(nx1:nx2, ny1:ny2) = sig*A;
            SIGzz(SIGzz <= 1)       = 1;

            if plotControl
                % SHOW DEVICE ERzz
                subplot(1, 2, 1);
                imagesc(xa, ya, ERzz.');
                axis equal tight off;
                title('ERzz');
                colorbar;
                plot_darkmode

                % SHOW DEVICE SIGzz
                subplot(1, 2, 2);
                imagesc(xa, ya, SIGzz.');
                axis equal tight off;
                title('SIGzz');
                colorbar;
                plot_darkmode;

                drawnow;
            end 

            indexName                    = (rot_dim + ldim1*step_var_pent - step_var_pent);

            imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
            imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
            Shape13(indexName).DevER      = imageER_zz ;
            Shape13(indexName).DevSIG     = imageSIG_zz ;
            Shape13(indexName).Label      = classNum(13);

            indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                            step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                            step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                            indexName;

            Shape13(indexName).Image_num  = indexLabel;

            disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(13)) ' Label ' num2str(indexLabel)])
        end
    end
    rng default
    ShapeCollection = [Shape13];
    Centered_Pentagon = [Shape13];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Centered Hexagon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if centerHex
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Random Values for R
    
    limitA = L/10;
    limitB = L/2;
    rng(seedControl_hex);
    RAD1x = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD2x = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD3x = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD4x = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD5x = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD6x = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD1y = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD2y = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD3y = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD4y = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD5y = limitA + (limitB - limitA).*rand(1,step_var_hex);
    RAD6y = limitA + (limitB - limitA).*rand(1,step_var_hex);

    
    % Handle Values
    rot_crs = linspace(0,90,step_var_hex)*degrees;

    for ldim1 = 1:length(RAD1x)             
        for rot_dim = 1:length(rot_crs)
            [TY, TX] = meshgrid(ya1, xa1);
            % ROTATE COORDINATES

            phi1            = rot_crs(rot_dim);   % handles
            [THETA, RADIUS] = cart2pol(TX, TY);
            [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

            % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

            p1 = [RAD1x(ldim1)*cos(phi1)              ; RAD1y(ldim1)*sin(phi1)];
            p2 = [RAD2x(ldim1)*cos(phi1+60*degrees)   ; RAD2y(ldim1)*sin(phi1+60*degrees)];
            p3 = [RAD3x(ldim1)*cos(phi1+60*2*degrees) ; RAD3y(ldim1)*sin(phi1+60*2*degrees)];
            p4 = [RAD4x(ldim1)*cos(phi1+60*3*degrees) ; RAD4y(ldim1)*sin(phi1+60*3*degrees)];
            p5 = [RAD5x(ldim1)*cos(phi1+60*4*degrees) ; RAD5y(ldim1)*sin(phi1+60*4*degrees)];
            p6 = [RAD6x(ldim1)*cos(phi1+60*5*degrees) ; RAD6y(ldim1)*sin(phi1+60*5*degrees)];
            
            P  = [p1 p2 p3 p4 p5 p6];

            A = polyfill(xa1,ya1,P);
            A = fliplr(A); % use fliplr rather tha flipud

            % STOP/START INDICES
            nx = round(L/dx);
            nx1 = 1 + floor((Nx - nx)/2);
            nx2 = nx1 + nx - 1;

            ny = round(L/dy);
            ny1 = 1 + floor((Ny - ny)/2);
            ny2 = ny1 + ny - 1;

            % INCORPORATE MATERIAL PROPERTIES
            ERzz(nx1:nx2, ny1:ny2)  = er*A;
            ERzz(ERzz <= 1)         = 1;
            SIGzz(nx1:nx2, ny1:ny2) = sig*A;
            SIGzz(SIGzz <= 1)       = 1;

            if plotControl
                % SHOW DEVICE ERzz
                subplot(1, 2, 1);
                imagesc(xa, ya, ERzz.');
                axis equal tight off;
                title('ERzz');
                colorbar;
                plot_darkmode

                % SHOW DEVICE SIGzz
                subplot(1, 2, 2);
                imagesc(xa, ya, SIGzz.');
                axis equal tight off;
                title('SIGzz');
                colorbar;
                plot_darkmode;

                drawnow;
            end 

            indexName                    = (rot_dim + ldim1*step_var_hex - step_var_hex);

            imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
            imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
            Shape14(indexName).DevER      = imageER_zz ;
            Shape14(indexName).DevSIG     = imageSIG_zz ;
            Shape14(indexName).Label      = classNum(14);

            indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                            step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                            step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                            step_var_pent^2 + indexName;

            Shape14(indexName).Image_num  = indexLabel;

            disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(14)) ' Label ' num2str(indexLabel)])
        end
    end
    rng default
    ShapeCollection = [Shape14];
    Centered_Hexagon = [Shape14];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Parallelogram Ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parallelogram_ring
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
   
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Random Values for R
    
    limitA = L/10;
    limitB = L/2;
    rng(seedControl_par_ring);
    
    RAD1x = limitA + (limitB - limitA).*rand(1,step_par_ring);
    RAD2x = limitA + (limitB - limitA).*rand(1,step_par_ring);
    RAD3x = limitA + (limitB - limitA).*rand(1,step_par_ring);
    RAD4x = limitA + (limitB - limitA).*rand(1,step_par_ring);
    
    RAD1y = limitA + (limitB - limitA).*rand(1,step_par_ring);
    RAD2y = limitA + (limitB - limitA).*rand(1,step_par_ring);
    RAD3y = limitA + (limitB - limitA).*rand(1,step_par_ring);
    RAD4y = limitA + (limitB - limitA).*rand(1,step_par_ring);
    
    % Handle Values
    rot_crs = linspace(0,180,step_par_ring)*degrees;
    scaling = linspace(1.5,5,step_par_ring);
    
    for scl = 1:length(scaling)   
        for ldim1 = 1:length(RAD1x)             
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);
                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

                p1 = [RAD1x(ldim1)*cos(phi1)              ; RAD1y(ldim1)*sin(phi1)];
                p2 = [RAD2x(ldim1)*cos(phi1+90*degrees)   ; RAD2y(ldim1)*sin(phi1+90*degrees)];
                p3 = [RAD3x(ldim1)*cos(phi1+90*2*degrees) ; RAD3y(ldim1)*sin(phi1+90*2*degrees)];
                p4 = [RAD4x(ldim1)*cos(phi1+90*3*degrees) ; RAD4y(ldim1)*sin(phi1+90*3*degrees)];

                P1  = [p1 p2 p3 p4];
                P2  = [p1 p2 p3 p4]/scaling(scl);

                A1 = polyfill(xa1,ya1,P1);
                A1 = fliplr(A1); % use fliplr rather tha flipud

                A2 = polyfill(xa1,ya1,P2);
                A2 = fliplr(A2); % use fliplr rather tha flipud

                A = double(xor(A1,A2));

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*A;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*A;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim1*step_par_ring - step_par_ring) + scl*step_par_ring^2 - step_par_ring^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape15(indexName).DevER      = imageER_zz ;
                Shape15(indexName).DevSIG     = imageSIG_zz ;
                Shape15(indexName).Label      = classNum(15);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                                step_var_pent^2 +  step_var_hex^2 + indexName;

                Shape15(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(15)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    rng default
    ShapeCollection = [Shape15];
    Parallelogram_Ring = [Shape15];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Triangle Ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if triangle_ring
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Random Values for R
    
    limitA = L/10;
    limitB = L/2;
    rng(seedControl_tri_ring);
    
    RAD1x = limitA + (limitB - limitA).*rand(1,step_tri_ring);
    RAD2x = limitA + (limitB - limitA).*rand(1,step_tri_ring);
    RAD3x = limitA + (limitB - limitA).*rand(1,step_tri_ring);
   
    RAD1y = limitA + (limitB - limitA).*rand(1,step_tri_ring);
    RAD2y = limitA + (limitB - limitA).*rand(1,step_tri_ring);
    RAD3y = limitA + (limitB - limitA).*rand(1,step_tri_ring);
   
    
    % Handle Values
    rot_crs = linspace(0,180,step_tri_ring)*degrees;
    scaling = linspace(1.5,5,step_tri_ring);
    
    for scl = 1:length(scaling)   
        for ldim1 = 1:length(RAD1x)             
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);
                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

                p1 = [RAD1x(ldim1)*cos(phi1)               ; RAD1y(ldim1)*sin(phi1)];
                p2 = [RAD2x(ldim1)*cos(phi1+120*degrees)   ; RAD2y(ldim1)*sin(phi1+120*degrees)];
                p3 = [RAD3x(ldim1)*cos(phi1+120*2*degrees) ; RAD3y(ldim1)*sin(phi1+120*2*degrees)];
               
                P1  = [p1 p2 p3];
                P2  = [p1 p2 p3]/scaling(scl);

                A1 = polyfill(xa1,ya1,P1);
                A1 = fliplr(A1); % use fliplr rather tha flipud

                A2 = polyfill(xa1,ya1,P2);
                A2 = fliplr(A2); % use fliplr rather tha flipud

                A = double(xor(A1,A2));

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*A;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*A;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim1*step_tri_ring - step_tri_ring) + scl*step_tri_ring^2 - step_tri_ring^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape16(indexName).DevER      = imageER_zz ;
                Shape16(indexName).DevSIG     = imageSIG_zz ;
                Shape16(indexName).Label      = classNum(16);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                                step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + indexName;

                Shape16(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(16)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    rng default
    ShapeCollection = [Shape16];
    Triangle_Ring = [Shape16];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Parallelogram with Circular Hole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parallelogram_hole
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Random Values for R
    
    limitA = L/6;
    limitB = L/2;
    rng(seedControl_par_hole);
    
    RAD1x = limitA + (limitB - limitA).*rand(1,step_par_hole);
    RAD2x = limitA + (limitB - limitA).*rand(1,step_par_hole);
    RAD3x = limitA + (limitB - limitA).*rand(1,step_par_hole);
    RAD4x = limitA + (limitB - limitA).*rand(1,step_par_hole);
    
    RAD1y = limitA + (limitB - limitA).*rand(1,step_par_hole);
    RAD2y = limitA + (limitB - limitA).*rand(1,step_par_hole);
    RAD3y = limitA + (limitB - limitA).*rand(1,step_par_hole);
    RAD4y = limitA + (limitB - limitA).*rand(1,step_par_hole);
    
    % Handle Values
    rot_crs  = linspace(0,180,step_par_hole)*degrees;
    hole_rad = linspace(L/10,L/7,step_par_hole);
    
    for rh = 1:length(hole_rad)   
        for ldim1 = 1:length(RAD1x)             
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);
                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

                p1 = [RAD1x(ldim1)*cos(phi1)              ; RAD1y(ldim1)*sin(phi1)];
                p2 = [RAD2x(ldim1)*cos(phi1+90*degrees)   ; RAD2y(ldim1)*sin(phi1+90*degrees)];
                p3 = [RAD3x(ldim1)*cos(phi1+90*2*degrees) ; RAD3y(ldim1)*sin(phi1+90*2*degrees)];
                p4 = [RAD4x(ldim1)*cos(phi1+90*3*degrees) ; RAD4y(ldim1)*sin(phi1+90*3*degrees)];

                P1  = [p1 p2 p3 p4];
                P2  = (TX.^2 + TY.^2) <= hole_rad(rh).^2;

                A1 = polyfill(xa1,ya1,P1);
                A1 = fliplr(A1); % use fliplr rather tha flipud

                A = double(xor(P2,A1));

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*A;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*A;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim1*step_par_hole - step_par_hole) + rh*step_par_hole^2 - step_par_hole^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape17(indexName).DevER      = imageER_zz ;
                Shape17(indexName).DevSIG     = imageSIG_zz ;
                Shape17(indexName).Label      = classNum(17);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                                step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                                indexName;

                Shape17(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(17)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    rng default
    ShapeCollection = [Shape17];
    Parallelogram_with_Hole = [Shape17];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Triangle with Circular Hole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if triangle_hole
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Random Values for R
    
    limitA = L/4;
    limitB = L/2;
    rng(seedControl_tri_hole);
    
    RAD1x = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    RAD2x = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    RAD3x = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    RAD4x = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    
    RAD1y = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    RAD2y = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    RAD3y = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    RAD4y = limitA + (limitB - limitA).*rand(1,step_tri_hole);
    
    % Handle Values
    rot_crs  = linspace(0,180,step_tri_hole)*degrees;
    hole_rad = linspace(L/12,L/9,step_tri_hole);
    
    for rh = 1:length(hole_rad)   
        for ldim1 = 1:length(RAD1x)             
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);
                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

                p1 = [RAD1x(ldim1)*cos(phi1)              ; RAD1y(ldim1)*sin(phi1)];
                p2 = [RAD2x(ldim1)*cos(phi1+120*degrees)   ; RAD2y(ldim1)*sin(phi1+120*degrees)];
                p3 = [RAD3x(ldim1)*cos(phi1+120*2*degrees) ; RAD3y(ldim1)*sin(phi1+120*2*degrees)];
                

                P1  = [p1 p2 p3];
                P2  = (TX.^2 + TY.^2) <= hole_rad(rh).^2;

                A1 = polyfill(xa1,ya1,P1);
                A1 = fliplr(A1); % use fliplr rather tha flipud

                A = double(xor(P2,A1));

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*A;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*A;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim1*step_tri_hole - step_tri_hole) + rh*step_tri_hole^2 - step_tri_hole^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape18(indexName).DevER      = imageER_zz ;
                Shape18(indexName).DevSIG     = imageSIG_zz ;
                Shape18(indexName).Label      = classNum(18);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                                step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                                step_par_hole^3 + indexName;

                Shape18(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(18)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    rng default
    ShapeCollection = [Shape18];
    Triangle_with_Hole = [Shape18];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Polygons - Regular Hexagon with Circular Hole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if reg_hex_hole
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Handle Values
    rot_crs   = linspace(0,180,step_hex_hole)*degrees;
    spoke_rad = linspace(L/8,L/2,step_hex_hole);
    hole_rad  = spoke_rad/2;
    
    for rh = 1:length(hole_rad)   
        for ldim1 = 1:length(spoke_rad)             
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);
                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

                p1 = [spoke_rad(ldim1)*cos(phi1)               ; spoke_rad(ldim1)*sin(phi1)];
                p2 = [spoke_rad(ldim1)*cos(phi1+60*degrees)    ; spoke_rad(ldim1)*sin(phi1+60*degrees)];
                p3 = [spoke_rad(ldim1)*cos(phi1+60*2*degrees)  ; spoke_rad(ldim1)*sin(phi1+60*2*degrees)];
                p4 = [spoke_rad(ldim1)*cos(phi1+60*3*degrees)  ; spoke_rad(ldim1)*sin(phi1+60*3*degrees)];
                p5 = [spoke_rad(ldim1)*cos(phi1+60*4*degrees)  ; spoke_rad(ldim1)*sin(phi1+60*4*degrees)];
                p6 = [spoke_rad(ldim1)*cos(phi1+60*5*degrees)  ; spoke_rad(ldim1)*sin(phi1+60*5*degrees)];
               
                P1  = [p1 p2 p3 p4 p5 p6];
                P2  = (TX.^2 + TY.^2) <= hole_rad(rh).^2;

                A1 = polyfill(xa1,ya1,P1);
                A1 = fliplr(A1); % use fliplr rather tha flipud

                A = double(xor(P2,A1));

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*A;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*A;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + ldim1*step_hex_hole - step_hex_hole) + rh*step_hex_hole^2 - step_hex_hole^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape19(indexName).DevER      = imageER_zz ;
                Shape19(indexName).DevSIG     = imageSIG_zz ;
                Shape19(indexName).Label      = classNum(19);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                                step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                                step_par_hole^3 + step_tri_hole^3 + indexName;

                Shape19(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(19)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape19];
    Reg_Hexagon_with_Hole = [Shape19];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Moon-shaped Objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if moon_shaped
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Handle Values
    rot_crs  = linspace(0,180,step_moon)*degrees;
    radius1  = linspace(L/4,L/2,step_moon);
    radius2  = linspace(L/8,L/3,step_moon);
    
    for rad2 = 1:length(radius2)   
        for rad1 = 1:length(radius1)             
            for rot_dim = 1:length(rot_crs)
                [TY, TX] = meshgrid(ya1, xa1);
                % ROTATE COORDINATES

                phi1            = rot_crs(rot_dim);   % handles
                [THETA, RADIUS] = cart2pol(TX, TY);
                [TX, TY]        = pol2cart(THETA+phi1, RADIUS);

                % PRECISE VERTICES (Translate Cartesian to Polar to allow Rotation)

                P1  = (TX.^2 + TY.^2) <= radius1(rad1).^2;
                P2  = ((TX - radius1(rad1)/1.5).^2 +...
                      (TY - radius1(rad1)/1.5).^2) <= radius2(rad2).^2;

                A = double(xor(P1,P1 & P2));

                % STOP/START INDICES
                nx = round(L/dx);
                nx1 = 1 + floor((Nx - nx)/2);
                nx2 = nx1 + nx - 1;

                ny = round(L/dy);
                ny1 = 1 + floor((Ny - ny)/2);
                ny2 = ny1 + ny - 1;

                % INCORPORATE MATERIAL PROPERTIES
                ERzz(nx1:nx2, ny1:ny2)  = er*A;
                ERzz(ERzz <= 1)         = 1;
                SIGzz(nx1:nx2, ny1:ny2) = sig*A;
                SIGzz(SIGzz <= 1)       = 1;

                if plotControl
                    % SHOW DEVICE ERzz
                    subplot(1, 2, 1);
                    imagesc(xa, ya, ERzz.');
                    axis equal tight off;
                    title('ERzz');
                    colorbar;
                    plot_darkmode

                    % SHOW DEVICE SIGzz
                    subplot(1, 2, 2);
                    imagesc(xa, ya, SIGzz.');
                    axis equal tight off;
                    title('SIGzz');
                    colorbar;
                    plot_darkmode;

                    drawnow;
                end 

                indexName                    = (rot_dim + rad1*step_moon - step_moon) + rad2*step_moon^2 - step_moon^2;

                imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
                imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
                Shape20(indexName).DevER      = imageER_zz ;
                Shape20(indexName).DevSIG     = imageSIG_zz ;
                Shape20(indexName).Label      = classNum(20);

                indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                                step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                                step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                                step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                                step_par_hole^3 + step_tri_hole^3 + step_hex_hole^3 + indexName;

                Shape20(indexName).Image_num  = indexLabel;

                disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(20)) ' Label ' num2str(indexLabel)])
            end
        end
    end
    ShapeCollection = [Shape20];
    Moon_shaped = [Shape20];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Hexagon (Square) Array of Circles with Center Circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hex_circ_array
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);

    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Handle Values
    rot_crs  = linspace(0,180,step_hex_circ)*degrees;
    radius1  = linspace(L/4,L/3,step_hex_circ);
      
    for rad1 = 1:length(radius1)             
        for rot_dim = 1:length(rot_crs)
            % ROTATE COORDINATES
            [Y2, X2] = meshgrid(ya1, xa1);

            phi1            = rot_crs(rot_dim);   % handles
            [THETA, RADIUS] = cart2pol(X2, Y2);
            [X2, Y2]        = pol2cart(THETA+phi1, RADIUS);

            a       = L;
            r       = radius1(rad1);
            NP      = 2;
            b       = L;

            % BUILD PHOTONIC CRYSTAL
            ER2 = zeros(Nx,Nx);
            y1  = -L/2;
            for np = 1 : NP + 1
                y0 = y1;
                ER2 = ER2 | ((X2 - a/2).^2 + (Y2 - y0).^2 < r^2);
                ER2 = ER2 | ((X2 + a/2).^2 + (Y2 - y0).^2 < r^2);
                ER2 = ER2 | ((X2 +   0).^2 + (Y2 - y0 - b/2).^2 < r^2);
                y1 = y1 + b;
            end

            %CLIP THE LATTICE
            y1  = NPML(1)*dy - SPACER(1);
            y2  = y1 + (b)*NP;
            %ER2 = 1 - ER2;
            ER2 = ER2.*(Y2 >= y1 & Y2 <=y2);
            ER2 = double(ER2);

            % STOP/START INDICES
            nx = round(L/dx);
            nx1 = 1 + floor((Nx - nx)/2);
            nx2 = nx1 + nx - 1;

            ny = round(L/dy);
            ny1 = 1 + floor((Ny - ny)/2);
            ny2 = ny1 + ny - 1;

            % INCORPORATE MATERIAL PROPERTIES
            ERzz(nx1:nx2, ny1:ny2)  = er*ER2;
            ERzz(ERzz <= 1)         = 1;
            SIGzz(nx1:nx2, ny1:ny2) = sig*ER2;
            SIGzz(SIGzz <= 1)       = 1;

            if plotControl
                % SHOW DEVICE ERzz
                subplot(1, 2, 1);
                imagesc(xa, ya, ERzz.');
                axis equal tight off;
                title('ERzz');
                colorbar;
                plot_darkmode

                % SHOW DEVICE SIGzz
                subplot(1, 2, 2);
                imagesc(xa, ya, SIGzz.');
                axis equal tight off;
                title('SIGzz');
                colorbar;
                plot_darkmode;

                drawnow;
            end 

            indexName                    = (rot_dim + rad1*step_hex_circ - step_hex_circ);

            imageER_zz                   = ERzz(1:nx_design,ny_1d:ny_2d);
            imageSIG_zz                  = SIGzz(1:nx_design,ny_1d:ny_2d);
            Shape21(indexName).DevER      = imageER_zz ;
            Shape21(indexName).DevSIG     = imageSIG_zz ;
            Shape21(indexName).Label      = classNum(21);

            indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                            step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                            step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                            step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                            step_par_hole^3 + step_tri_hole^3 + step_hex_hole^3 + step_moon^3 + ...
                                            indexName;

            Shape21(indexName).Image_num  = indexLabel;

            disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(21)) ' Label ' num2str(indexLabel)])
        end
    end
    ShapeCollection = [Shape21];
    Hex_Array_Circles = [Shape21];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Randomized Binary Design 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if binary_array
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Rows and Columns of the Binary Array
    e_w = num_elem_binary;
    e_h = e_w;
    
    % Initialize the single element space
    fill         = ones(ceil(Nx/e_w), ceil(Nx/e_h));
    
    % Initialize the Entire Device Space
    binary_array = zeros(Nx, Nx);
    
    % The Single Element Space
    [sing_row_el, sing_col_el]= size(fill);

    % Limits of the random values
    limitA = 0;
    limitB = 1;
    N = e_w.^2;
    
    for iter = 1:step_bin_array
      
        % Binary values
        bind_rand_values = randi([limitA, limitB], 1, N);

        c_row = 1;
        for row_index = 1:e_h
            c_col = 1;
            for col_index = 1:e_w

                binary_array(c_row:row_index*sing_row_el, c_col: col_index*sing_col_el) = fill*bind_rand_values(col_index + (e_w*(row_index - 1)));

                c_col = col_index*sing_col_el+1;
            end
            c_row = row_index*sing_row_el+1;
        end

        % STOP/START INDICES
        nx = round(L/dx);
        nx1 = 1 + floor((Nx - nx)/2);
        nx2 = nx1 + nx - 1;

        ny = round(L/dy);
        ny1 = 1 + floor((Ny - ny)/2);
        ny2 = ny1 + ny - 1;

        % INCORPORATE MATERIAL PROPERTIES
        ERzz(nx1:nx2, ny1:ny2)  = er*binary_array;
        ERzz(ERzz <= 1)         = 1;
        SIGzz(nx1:nx2, ny1:ny2) = sig*binary_array;
        SIGzz(SIGzz <= 1)       = 1;

        if plotControl
            % SHOW DEVICE ERzz
            subplot(1, 2, 1);
            imagesc(xa, ya, ERzz.');
            axis equal tight off;
            title('ERzz');
            colorbar;
            plot_darkmode

            % SHOW DEVICE SIGzz
            subplot(1, 2, 2);
            imagesc(xa, ya, SIGzz.');
            axis equal tight off;
            title('SIGzz');
            colorbar;
            plot_darkmode;

            drawnow;
        end 

        indexName                     = iter;
        imageER_zz                    = ERzz(1:nx_design,ny_1d:ny_2d);
        imageSIG_zz                   = SIGzz(1:nx_design,ny_1d:ny_2d);
        Shape22(indexName).DevER      = imageER_zz ;
        Shape22(indexName).DevSIG     = imageSIG_zz ;
        Shape22(indexName).Label      = classNum(22);

        indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                        step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                        step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                        step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                        step_par_hole^3 + step_tri_hole^3 + step_hex_hole^3 + step_moon^3 + ...
                                        step_hex_circ^2 + indexName;

        Shape22(indexName).Image_num  = indexLabel;

        disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(22)) ' Label ' num2str(indexLabel)])
    end
    ShapeCollection = [Shape22];
    Binary_Design = [Shape22];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE CONSTRUCTION: Randomized Pixel Design 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pixel_array
    
    ERzz    = ones(Nx,Ny);
    URxx    = ones(Nx,Ny);
    URyy    = ones(Nx,Ny);
    SIGzz   = zeros(Nx,Ny);
    
    % Design Space
    nx_design = round(L/dx);
    ny_design = round(L/dy);
    ny_1d = 1 + floor((Ny - ny_design)/2);
    ny_2d = ny_1d + ny_design - 1;
    
    % Create a Meshgrid of Nx-by-Nx dimension
    
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    % Rows and Columns of the Binary Array
    e_w = num_elem_cont;
    e_h = e_w;
    
    % Initialize the single element space
    fill         = ones(ceil(Nx/e_w), ceil(Nx/e_h));
    
    % Initialize the Entire Device Space
    pixel_array = zeros(Nx, Nx);
    
    % The Single Element Space
    [sing_row_el, sing_col_el]= size(fill);

    % Limits of the random values
    limitA = 0;
    limitB = 1;
    N = e_w.^2;
    
    for iter = 1:step_bin_array
        % Continuous values
        cont_rand_values = limitA + (limitB - limitA).*rand(1,N); 

        c_row = 1;
        for row_index = 1:e_h
            c_col = 1;
            for col_index = 1:e_w

                pixel_array(c_row:row_index*sing_row_el, c_col: col_index*sing_col_el) = fill*cont_rand_values(col_index + (e_w*(row_index - 1)));

                c_col = col_index*sing_col_el+1;
            end
            c_row = row_index*sing_row_el+1;
        end

        % STOP/START INDICES
        nx = round(L/dx);
        nx1 = 1 + floor((Nx - nx)/2);
        nx2 = nx1 + nx - 1;

        ny = round(L/dy);
        ny1 = 1 + floor((Ny - ny)/2);
        ny2 = ny1 + ny - 1;

        % INCORPORATE MATERIAL PROPERTIES
        ERzz(nx1:nx2, ny1:ny2)  = er*pixel_array;
        ERzz(ERzz <= 1)         = 1;
        SIGzz(nx1:nx2, ny1:ny2) = sig*pixel_array;
        SIGzz(SIGzz <= 1)       = 1;

        if plotControl
            % SHOW DEVICE ERzz
            subplot(1, 2, 1);
            imagesc(xa, ya, ERzz.');
            axis equal tight off;
            title('ERzz');
            colorbar;
            plot_darkmode

            % SHOW DEVICE SIGzz
            subplot(1, 2, 2);
            imagesc(xa, ya, SIGzz.');
            axis equal tight off;
            title('SIGzz');
            colorbar;
            plot_darkmode;

            drawnow;
        end 

        indexName                     = iter;
        imageER_zz                    = ERzz(1:nx_design,ny_1d:ny_2d);
        imageSIG_zz                   = SIGzz(1:nx_design,ny_1d:ny_2d);
        Shape23(indexName).DevER      = imageER_zz ;
        Shape23(indexName).DevSIG     = imageSIG_zz ;
        Shape23(indexName).Label      = classNum(23);

        indexLabel                    = step_var_rec^2 + step_var_circ + step_var_tri^3 + step_var_ell^2 +...
                                        step_var_crs^3 + step_var_dcrs^3 + step_var_rtri^2 + step_var_rrect^3 + ...
                                        step_var_rell^3 + step_var_comc^3 + step_var_ring^3 + step_var_trap^3 + ...
                                        step_var_pent^2 +  step_var_hex^2 + step_par_ring^3 + step_tri_ring^3 + ...
                                        step_par_hole^3 + step_tri_hole^3 + step_hex_hole^3 + step_moon^3 + ...
                                        step_hex_circ^2 + step_bin_array+ indexName;

        Shape23(indexName).Image_num  = indexLabel;

        disp(['Processing Image: ' num2str(indexName) ' Class: ' num2str(classNum(23)) ' Label ' num2str(indexLabel)])
    end
    
    ShapeCollection = [Shape23];
    Pixel_Design = [Shape23];
end            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE COLLECTION: CONCATENATE STRUCTURE ARRAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if collate
 
    ShapeCollection = [Shape1, Shape2, Shape3, Shape4, Shape5, Shape6, Shape7, ...
                       Shape8, Shape9, Shape10, Shape11, Shape12, Shape13, Shape14, ...
                       Shape15, Shape16, Shape17, Shape18, Shape19, Shape20, Shape21]; %Shape22, Shape23];
                   
    checkValue1     = sum(linspace(1, total_num_of_samples, total_num_of_samples));
    checkValue2     = sum([ShapeCollection(:).Image_num]);
    disp(['Size of the Shape Collection: ' '[' num2str(size(ShapeCollection)) ']'])
    
    if checkValue1 == checkValue2
        disp(['Check number of elements and numbering: ' '[' num2str(checkValue1) ' ' num2str(checkValue2) ']' ' ' 'Equal'])
    else
        disp(['Check number of elements and numbering: ' '[' num2str(checkValue1) ' ' num2str(checkValue2) ']' ' ' 'Not Equal'])
    end
    
    if save_the_shapes
        % SAVE THE SHAPES
        save(Filename,'ShapeCollection','-v7.3','-nocompression');
        disp('Saving the collection of shapes, DONE!')
    end
    
    if save_the_shapes_HDF5
        
        save_time_1 = clock;
        disp('Saving in HDF5 file format')
        [rows_mat, cols_mat] = size(Shape1(1).DevER);
        
        for index = 1:length(ShapeCollection)
            h5create(filename_hdf5,strcat('/material/DevER/', num2str(index-1)),[rows_mat cols_mat], 'ChunkSize',[20 20],'Deflate',9)
        end
        
        for file_index = 1:length(ShapeCollection)
            h5write(filename_hdf5,strcat('/material/DevER/', num2str(file_index-1)), ShapeCollection(file_index).DevER);
        end
        disp('Saving in HDF5 file format --- DONE')
        save_time_2 = clock;
        
        time_s = etime(save_time_2, save_time_1);
        disp(['HDF5 Saving Elapsed time is ' num2str(time_s) ' seconds.']);
        disp(['HDF5 Saving Elapsed time is ' num2str(time_s/60) ' minutes.']);
        disp(['HDF5 Saving Elapsed time is ' num2str(time_s/60/60) ' hours.']);
        disp(['HDF5 Saving Elapsed time is ' num2str(time_s/60/60/24) ' days.']);

    end
end

if pixel_binary_group
    pixel_binary_Collection = [Shape22, Shape23];
                   
    disp(['Size of the Shape Collection: ' '[' num2str(size(pixel_binary_Collection)) ']'])

    if save_the_shapes
        % SAVE THE SHAPES
        save(Filename_pixel_binary,'pixel_binary_Collection','-v7.3','-nocompression');
        disp('Saving the collection of shapes, DONE!')
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE COLLECTION: DISPLAY ENTIRE COLLECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if display_samples
    xa1 = [0:Nx-1]*dx;   xa1 = xa1 - mean(xa1);
    ya1 = [0:Nx-1]*dy;   ya1 = ya1 - mean(ya1);
    
    Array1 = [];
    example_height = example_width;
    if collate
        rand_samples = randi([1, size(ShapeCollection,2)], 1, example_width^2);
        for array_index = 1:(example_width^2)
           A = ShapeCollection(rand_samples(array_index)).DevER;
           Array1 = horzcat(Array1, A);
        end
    end
    
    if pixel_binary_group
        rand_samples = randi([1, size(pixel_binary_Collection,2)], 1, example_width^2);
        for array_index = 1:(example_width^2)
           A = pixel_binary_Collection(rand_samples(array_index)).DevER;
           Array1 = horzcat(Array1, A);
        end
    end
    
   
    [m n] = size(Array1);
    display_rows = floor(m * example_height);
    display_cols = floor(n / example_width);
    
    display_array = zeros(display_rows, display_cols);
    
    current_row = 1;
    for row_index = 1:example_height
        current_col = 1;
        for col_index = 1:example_width
            display_array(current_row:row_index*m, current_col:col_index*m) = ...
                    Array1(:, current_col + (display_cols*(row_index - 1)) ...
                    :col_index*m + (display_cols*(row_index - 1)));
            current_col = col_index*m+1;
        end
        current_row = row_index*m+1;
    end
    
    
    % SHOW DEVICE ERzz
    subplot(1, 2, 1);
    imagesc(xa1, ya1, display_array);
    axis equal tight off;
    title('ERzz');
    colorbar;
    plot_darkmode
   
    % 3D View ERzz
    subplot(1, 2, 2);
    mesh(display_array);
    title('ERzz 3D');
    colorbar;
    plot_darkmode
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE COLLECTION: SAVE INDIVIDUAL SHAPE TYPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if save_rec
    ShapeType = [Shape1]
    save('recShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_circ   
    ShapeType = [Shape2];
    save('circShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_tri
    ShapeType = [Shape3];
    save('triShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_ell   
    ShapeType = [Shape4];
    save('ellShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_crs
    ShapeType = [Shape5];
    save('crsShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_dcrs 
    ShapeType = [Shape6];
    save('dcrsShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_rtri
    ShapeType = [Shape7];
    save('rtriShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_rrect 
    ShapeType = [Shape8];
    save('rrectShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end
    
if save_rell 
    ShapeType = [Shape9];
    save('rellShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_comc  
    ShapeType = [Shape10];
    save('comcShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_ring  
    ShapeType = [Shape11];
    save('ringShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_trap 
    ShapeType = [Shape12];
    save('trapShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_pent 
    ShapeType = [Shape13];
    save('pentShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_hex   
    ShapeType = [Shape14];
    save('hexShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_par_ring 
    ShapeType = [Shape15];
    save('parRingShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_tri_ring  
    ShapeType = [Shape16];
    save('triRingShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_par_hole  
    ShapeType = [Shape17];
    save('parHoleShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_tri_hole  
    ShapeType = [Shape18];
    save('triHoleShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_hex_hole  
    ShapeType = [Shape19];
    save('hexHoleShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_moon  
    ShapeType = [Shape20];
    save('moonShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_hex_circ   
    ShapeType = [Shape21];
    save('hexCircShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_bin_array
    ShapeType = [Shape22];
    save('binaryShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

if save_pixel_arr  
    ShapeType = [Shape23];
    save('pixelShapes','ShapeType','-v7.3','-nocompression');
    disp('Saving the shapes, DONE!')
end

t2 = clock;
t = etime(t2,t1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);

