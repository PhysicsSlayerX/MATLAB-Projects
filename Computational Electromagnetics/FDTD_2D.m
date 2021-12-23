

numgraph_row = 1;
numgraph_col = 2;


er = 2.0;
ur = 3.0;

%Define Grid
Nx = 50;
Ny = 60;
dx = 1;
dy = 1;

% 2X Grid
Nx2 = 2*Nx;
Ny2 = 2*Ny;
dx2 = dx/2;
dy2 = dy/2;

%Create a Cylinder
r = 20;
xa2 = [0:Nx2-1]*dx2;
ya2 = [0:Ny2-1]*dy2;
xa2 = xa2 - mean(xa2);
ya2 = ya2 - mean(ya2);

[Y2, X2] = meshgrid(ya2, xa2);

ER = (X2.^2 + Y2.^2) <= r^2;
UR = (X2.^2 + Y2.^2) <= r^2;

ER2 = er.*ER;
UR2 = ur.*UR;
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

h = imagesc(xa2,ya2, ERxx.');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\epsilon_r');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;

subplot(numgraph_row,numgraph_col,2);

g = imagesc(xa2,ya2, URxx.');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('\mu_r');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;

drawnow; 


% Dirichlet Boundary Condition for CEx

for nx = 1:Nx
   for ny = 1 : Ny - 1
      CEx(nx, ny) = (Ez(nx, ny+1) - Ez(nx,ny))/dy; 
   end
   CEx(nx, Ny) = (0 - Ez(nx,ny))/dy;
end

% Dirichlet Boundary Condition for CEy

for ny = 1:Ny
   for nx = 1 : Nx - 1
      CEy(nx, ny) = -(Ez(nx+1, ny+1) - Ez(nx,ny))/dx; 
   end
   CEy(Nx, ny) = (0 - Ez(Nx,ny))/dx;
end

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
% 
