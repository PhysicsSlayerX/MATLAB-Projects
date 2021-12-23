close all
clc;
clear all;

fig = figure('Color','w');

degrees = pi/180;

k0 = 2*pi;
theta = 20*degrees;
kinc = k0*[sin(theta);cos(theta)];

dx = 1;
dy = 1;

xbc = -2;
ybc = -2; 

Nx = 3;
Ny = 3;

NGRID = [Nx Ny];
RES = [dx dy];
BC = [xbc ybc];

[DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc);
dex = full(DEX);
dey = full(DEY);
dhx = full(DHX);
dhy = full(DHY);

xa = [-Nx/2:Nx/2]*-dx;
ya = [0:Ny-1]*-dy;

subplot (1,4,1);
g = imagesc(xa,ya, real(dex).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DX');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));
subplot (1,4,2);
g = imagesc(xa,ya, real(dey).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DY');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));

subplot (1,4,3);
g = imagesc(xa,ya, real(dhx).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DX');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));
subplot (1,4,4);
g = imagesc(xa,ya, real(dhy).');
g2 = get(g, 'Parent');
set(g2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('GRID Derivative Matrix DY');
axis([-1 +1 +0 +10])
axis equal tight;
colorbar;
colormap(jet(1024));

set(gcf, 'Position',  [300, 200, 1500, 700])


%  disp(DEX);

%  disp(DHX);
%  disp(DHY);


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
        
        if BC(1) == 0
            DEX = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            diagonal_higher(1:NGRID(1):NGRID(1)*NGRID(2)) = 0;
            DEX = spdiags([-diagonal_0th diagonal_higher], [0 1], DEX);
            DEX = DEX/RES(1);
        elseif BC(1) == -2
            DEX = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            diagonal_higher(1:NGRID(1):NGRID(1)*NGRID(2)) = 0;
            diagonal_lower = zeros(n,1);
            periodicityX = NGRID(1)*RES(1);
            diagonal_lower(1:NGRID(1):NGRID(1)*NGRID(2)) = exp(1i*kinc(1)*periodicityX);  
            DEX = spdiags([diagonal_lower -diagonal_0th diagonal_higher], [-NGRID(1)+1 0 1], DEX);
            DEX = DEX/RES(1);
        end  
        
    end

    if NGRID(2) == 1
        I = eye();
        DEY = 1i*kinc(2)*I;

    else
       n = NGRID(1)*NGRID(2);
       
       if BC(2) == 0
            DEY = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            DEY = spdiags([-diagonal_0th diagonal_higher], [0 NGRID(1)], DEY);
            DEY = DEY/RES(2);
        elseif BC(2) == -2
            DEY = sparse(n,n);
            diagonal_0th = ones(n,1);
            diagonal_higher = ones(n,1);
            diagonal_lower = zeros(n,1);
            periodicityY = NGRID(2)*RES(2);
            diagonal_lower(1:1:NGRID(1)*NGRID(2)) = exp(1i*kinc(2)*periodicityY);  
            DEY = spdiags([diagonal_lower -diagonal_0th diagonal_higher], [-(NGRID(1)*NGRID(2)-(NGRID(1)-1))+1 0 NGRID(1)], DEY);
            DEY = DEY/RES(2);
       end  
       
    end
    
    DHX = -(DEX)';
    DHY = -(DEY)';
    
end