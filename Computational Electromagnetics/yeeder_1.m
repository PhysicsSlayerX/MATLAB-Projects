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
    
    if NGRID(1) > NGRID(2)
        indexX = NGRID(2);
        indexY = NGRID(1);
        
    elseif NGRID(1) == NGRID(2)
        indexX = NGRID(2);
        indexY = NGRID(1);
    else
        indexX = NGRID(2);
        indexY = NGRID(1);
    end
    
    if NGRID(1) == 1
        I = eye();
        DEX = 1i*kinc(1)*I;
        
    else
        n = NGRID(1)*NGRID(2);
        value = ones(n,1);
        DEX = (1/RES(1))*spdiags([-value value], [0 1], n,n);
        for a = 1:indexX
            if BC(1) == 0 && NGRID(1) < NGRID(2)
                
                DEX(a*indexX, a*indexX+1) = 0;
                
                if a >= indexX - 1
                    break;
                end
                
            elseif BC(1) == 0 && NGRID(1) > NGRID(2)
                DEX(a*NGRID(2), a*NGRID(2)+1) = 0;
                
            elseif BC(1) == 0 && NGRID(1) == NGRID(2)
                DEX(a*indexX, a*indexX+1) = 0;
                
            elseif BC(1) == -2
                DEX(a*NGRID(1), a*NGRID(1)+1) = 0;
                periodicityX = NGRID(1)*RES(1);
                DEX(a*NGRID(1), NGRID(1)*(a-1)+1) = exp(1i*kinc(1)*periodicityX);

            end
        end
    end

    if NGRID(2) == 1
        I = eye();
        DEY = 1i*kinc(2)*I;

    else
        n = NGRID(1)*NGRID(2);
        value = ones(n,1);
        DEY = (1/RES(2))*spdiags([-value value], [0 NGRID(1)], n,n);
        for a = 1:indexY
            if BC(2) == 0
                break;
                %DEY(a*indexY, a*indexY+1) = 0;
            elseif BC(2) == -2
                periodicityY = NGRID(2)*RES(2);
                DEY(NGRID(1)*(NGRID(2)-1)+a, a) = exp(1i*kinc(2)*periodicityY);                
            end
        end
    end
    if BC(1) == 0 && NGRID(1) > NGRID(2)
        DHX = -(DEX)';
        DHY = -(DEY)';
        
    else
        DEX(:,NGRID(1)*NGRID(2)+1) = [];
        DHX = -(DEX)';
        DHY = -(DEY)';
    end
    
end




