



format short

Nx = 7;
Ny = 4;
Nxlo = 2;
Nxhi = 3;
Nylo = 1;
Nyhi = 2;

NGRID = [Nx Ny];
NPML = [Nxlo Nxhi Nylo Nyhi];


[sx, sy] = calcpml2d(NGRID,NPML);

disp(sx);
disp(sy);

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
    
    a_max = 5;
    p = 5;
    sigmaprime_max = 1;
    eta0 = 376.73;

    sx = ones(NGRID);
    sy = sx ;
    
    for nx = 1:NPML(1)
        sx(NPML(1)-nx+1,:) = (1 + a_max*(nx/NPML(1))^p)*(1 + 1i*eta0*sigmaprime_max*(sin((pi*nx)/(2*NPML(1))))^2);

    end

    for nx = 1 : NPML(2)
        sx(NGRID(1)-NPML(2)+nx, :) = (1 + a_max*(nx/NPML(2))^p)*(1 + 1i*eta0*sigmaprime_max*(sin((pi*nx)/(2*NPML(2))))^2);

    end

    for ny = 1:NPML(3)
        sy(:,NPML(3)-ny+1) = (1 + a_max*(ny/NPML(3))^p)*(1 + 1i*eta0*sigmaprime_max*(sin((pi*ny)/(2*NPML(3))))^2);

    end

    for ny = 1 : NPML(4)
        sy(:,NGRID(2)-NPML(3)+ny) = (1 + a_max*(ny/NPML(3))^p)*(1 + 1i*eta0*sigmaprime_max*(sin((pi*ny)/(2*NPML(3))))^2);

    end
end