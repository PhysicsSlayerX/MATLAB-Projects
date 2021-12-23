%%Visualization of Matrix Operations

close all;
clc;
clear all;

%Matrices
A = [1 2 3; 4 5 6; 7 8 9; 10 11 12]

%Vector
v = [1;2;3]

%Size of the vector
[m n] = size(A)

%Another storing method
dim_A = size(A)

%Dimension of vector v
dim_v = size(v)

%Now let's index into the 2nd row 3rd column of Matrix A
A_23 = A(2,3)