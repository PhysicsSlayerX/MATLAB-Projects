% DEFINE GRID

Sx = 1;
Sy = 1;
Nx = 200;
Ny = 200;

% GRID ARRAYS

dx = Sx/Nx;
xa = [0:Nx - 1]*dx;
xa = xa - mean(xa);

dy = Sy/Ny;
ya = [0:Ny - 1]*dy;
ya = ya - mean(ya);

[Y, X] = meshgrid(ya, xa);

% DEFINE TWO POINTS
p1 = [0.50, 0.00];
p2 = [-0.50, -0.50];
p1 = [0.50, 0.00];
p2 = [-0.50, -0.50];

% FILL HALF SPACE
m = (p2(2)-p1(2))/(p2(1) - p1(1));
A = (Y - p1(2)) - m*(X  - p1(1)) > 0;
A = fliplr(A);

% PLOT A
subplot(1, 1, 1);
imagesc(xa, ya, A.');
axis equal tight;
title('ERzz');
colorbar;
plot_darkmode
drawnow;