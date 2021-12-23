clc;
clear;

Sx = 1;
Sy = 1;
Nx = 20;
Ny = 20;

dx = Sx/Nx;
xa = [0:Nx-1]*dx;
xa = xa - mean(xa);

dy = Sy/Ny;
ya = [0:Ny-1]*dy;
ya = ya - mean(ya);

[Y, X] = meshgrid(ya,xa);

x1 = -0.50;
y1 = +0.25;
x2 = +0.50;
y2 = -0.25;

m = (y2 - y1)/(x2 - x1);
A = (Y - y1) -  m*(X - x1) > 0;

fig = figure('Color','w');
set(gcf, 'Position', [1921 41 1920 963]);

h = imagesc(xa,ya, A.'); %shows the cartesian plane in the graph but not the ticks
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\mu_r');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;