%demo_graphics2d

%Initialize MATLAB

close all;
clc;
clear all;

% rand('seed' , 0);
% A = rand(200,200)
% 
% pcolor(A);
% shading interp;
% axis equal tight
% colorbar;
% 

% FUNCTION SIZE

Sx = 2;
Sy = 2;

Nx = 128;
Ny = round(Nx*Sy/Sx);

% GENERATE RANDOM FUNCTION
f = 0.02;
rand('seed', 0);

nx1 = 1 + floor(0.9999*f*Nx);
nx2 = Nx - nx1 + 1;

ny1 = 1 + floor(0.9999*f*Ny);
ny2 = Ny - ny1 + 1;

F               = rand(Nx,Ny) - 0.5;
F               = fft2(F);
F(nx1:nx2, :)   = 0;
F(:,ny1:ny2)    = 0; 
F               = real(ifft2(F));

% FUNCTION AXES

dx = Sx/Nx;
xa = [0.5:Nx-0.5]*dx;
dy = Sy/Ny;
ya = [0.5:Ny-0.5]*dy;

pcolor(xa, ya, F.');
shading interp;
set(gca, 'YDir' , 'normal');

% SET VIEW

axis equal tight;
set(gca, 'FontSize', 12, 'LineWidth', 4);
title('Weird Function $f(x,y)$', ...
    'Interpreter', 'Latex',  'FontSize', 22);

% X AXIS

xlim([0 Sx]);
T = [0: 0.5:Sx];
L = {};
for m = 1 : length(T)
    if T(m) == 0
        L{m} = 0;
    else
    L{m} = num2str(T(m), '%3.1f');
    end
end
set(gca, 'XTick', T, 'XTickLabel', L);
xlabel('$x$-axis','Interpreter', 'LaTeX', 'FontSize', 22);


ylim([0 Sy]);
T = [0: 0.5:Sy];
L = {};
for m = 1 : length(T)
    if T(m) == 0
        L{m} = 0;
    else
    L{m} = num2str(T(m), '%3.1f');
    end
end
set(gca, 'YTick', T, 'YTickLabel', L);
ylabel('$y$-axis','Interpreter', 'LaTeX', 'FontSize', 22);


