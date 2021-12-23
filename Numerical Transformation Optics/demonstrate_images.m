% demonstrate_images.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%imshow(images{1});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE A GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID PARAMETERS
Sx = 10;
Sy = 10;

Nx = 100;
Ny = round(Nx*Sy/Sx);

% CALCULATE GRID RESOLUTION
dx = Sx/Nx;
dy = Sy/Ny;

% CALCULATE AXIS VECTORS
xa = [0.5:Nx-0.5]*dx;
ya = [0.5:Ny-0.5]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD IMAGE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list of all PNG files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
addpath 'D:\Research\Matlab practice\Numerical Transformation Optics\circles'
imagefiles = dir('D:\Research\Matlab practice\Numerical Transformation Optics\circles\*.png');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   images{ii} = currentimage;
end
thisImage = imresize(images{end},[Nx Ny]);
lastimage = ceil(thisImage(:,:,1)/255);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY ARRAY A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW A
imagesc(xa, ya, lastimage.');

% SET VIEW
axis equal tight;
colorbar;





















