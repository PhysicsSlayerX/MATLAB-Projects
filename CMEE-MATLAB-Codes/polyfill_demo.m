% polyfill_demo.m
%
% POLYFILL      Fill 2D Grid with a Polygon
% 
% A = polyfill(xa,ya,P);
%
% xa,ya     Grid Axes for the Array A
% P         List of vertices for the polygon
%           [ x1 x2 ... xN ;
%             y1 y2 ... yN ];   N vertices
% A         2D array with polygon filled
%
% Note: The list of points P should progress CCW around the polygon.
%       If the points are listed in reverse order, the outside of
%       the polygon will be filled.

% INITIALIZE MATLAB
close all;
clc;
clear all;

% GRID
Sx = 1;
Sy = 1;
Nx = 64;
Ny = round(Nx*Sy/Sx);
xa = linspace(0,Sx,Nx);
ya = linspace(0,Sy,Ny);

% CREATE ARBITRARY POLYGON
p1 = [ 0.3 ; 0.1 ];
p2 = [ 0.8 ; 0.2 ];
p3 = [ 0.7 ; 0.9 ];
p4 = [ 0.6 ; 0.4 ];
p5 = [ 0.1 ; 0.8 ];
P  = [ p1 p2 p3 p4 p5 p1 ];

% CALL POLYFILL TO FILL POLYGON IN ARRAY A
A = polyfill(xa,ya,P);

% SHOW A
imagesc(xa,ya,A');
axis equal tight;
colorbar;












