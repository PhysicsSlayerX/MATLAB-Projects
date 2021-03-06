function [DX, D2X, DY, D2Y] = toder(NS, RES)
% TODER     Transformation Optics Derivative Matrices
%
% [DX, D2X, DY, D2Y] = toder(NS, RES);
%
% This MATLAB dunction generates derivative matrices for use in numerical
% transformation optics. Dirchlet boundary conditions are used.
%
% INPUT ARGUMENTS
% ===============
% NS    [Nx Ny] Size of the Grid
% RES   [dx dy] grid resolution parameters
%
% OUTPUT ARGUMENTS
% ===============
% DX    1st-order derivative matrix with respect to x
% D2X   2nd-order derivative matrix with respect to x
% DY    1st-order derivative matrix with respect to y
% D2Y   2nd-order derivative matrix with respect to y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT AND OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT SIZE OF THE GRID
Nx = NS(1);
Ny = NS(2);

% COMPUTEE SIZE OF MATRICES
M = Nx*Ny;

% COMPUTE SPECIAL MATRICES
Z = sparse(M, M);
V = ones(M, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DX and D2X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATRICES
DX = Z;
D2X = Z;

% BUILD DX and D2X
if Nx > 1
    
    % Build Default DX
    DX = spdiags(-V,-1,DX);
    DX = spdiags(V,+1,DX);
    
    % Build Default D2X
    D2X = spdiags(  +V,-1,D2X);
    D2X = spdiags(-2*V, 0,D2X);
    D2X = spdiags(  +V,+1,D2X);
    
    % Fix Boundary Erros (Default to Dirichlet)
    for ny = 1:Ny
        mxlo = (ny - 1)*Nx + 1;
        mxhi = (ny - 1)*Nx + Nx;
        if mxlo > 1
            DX(mxlo, mxlo-1) = 0;
            D2X(mxlo, mxlo-1) = 0;
        end
        if mxhi < M
            DX(mxhi, mxhi+1) = 0;
            D2X(mxhi, mxhi+1) = 0;
        end
    end
    
    % Divide by Grid Resolution
    DX = DX / (2 * RES(1));
    D2X = D2X / RES(1)^2;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DY and D2Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATRICES
DY = Z;
D2Y = Z;

% BUILD DX and D2X
if Ny > 1
    
    % Build DY
    DY = spdiags(-V,-Nx,DY);
    DY = spdiags(+V,+Nx,DY);
    
    % Build DY
    D2Y = spdiags(  +V, -Nx,D2Y);
    D2Y = spdiags(-2*V,   0,D2Y);
    D2Y = spdiags(  +V, +Nx,D2Y);
    
    % Divide by Grid Resolution
    DY = DY / (2 * RES(2));
    D2Y = D2Y / RES(2)^2;
    
end










































