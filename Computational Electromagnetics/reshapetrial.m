
clc;
clear;

Nx = 5;
Ny = 5;

Q = zeros(Nx, Ny);
n = Nx*Ny;
sr = 1;
Q(:, 1:2) = sr
Q1 = reshape(Q.',[],1)
Q2 = spdiags([Q1], [ 0 ], n,n);