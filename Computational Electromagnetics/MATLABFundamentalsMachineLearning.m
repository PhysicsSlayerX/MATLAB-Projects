%% MATLAB Fundamentals for Andrew Ng's Machine Learning Course

close all;
clc;
clear all;
a = pi;
disp(sprintf('2 decimals: %0.2f', a))

disp(sprintf('6 decimals: %0.6f', a))

OneMatrix = ones(2, 3)

ZeroMatrix = zeros(1,3)

randomMatrix = rand(3,3)

gaussianRand = randn(5,5)

identityMatrix = eye(6)

w = -6 + sqrt(10)*(randn(1, 10000));
%hist(w,50);

%% Moving Data Around 

A = [1 2; 3 4; 5 6]
sizeA = size(A)
sizeOfsize = size(sizeA)
sizeRow = size(A,1)
sizeColumn = size(A,2)
V = [1 2 3 4]
lengthV = length(V)
lengthA = length(A)
lengthSomeVector = length([1;2;3;4;5])

% pwd %current path or folder
% cd D:\ %change to D drive
% pwd %we're in the D drive
% ls %tells the contents of the folder/directory
% cd ('D:\Research\Matlab practice\Computational Electromagnetics')
% ls


T = readtable('FSSParametricData.csv','PreserveVariableNames',true);
Data = table2array(T);
GapData = Data(:,1);
RotationData = Data(:,2);

size(GapData);
size(RotationData);

%whos

%save GapData.mat GapData; %save the GapData vector into GapData.mat file

MatrixA = [1 2; 3 4; 5 6];
A_32 = A(3,2)
A_row2elements = A(2,:)
A_column2elements = A(:,2)

B = A([1 3],:)

A(:,2) = [10; 11; 12]

A = [A, [100;101;102]] %Adds another column vector

A(:) %put all elements of A into a single vector

A = [1 2; 3 4; 5 6]
B = [11 12;13 14; 15 16]

C = [A B] %concatenating two matrices (side by side)

D = [A; B] %concatenating two matrices (top bottom)

%%Computational

A = [1 2; 3 4; 5 6]
B = [11 12; 13 14; 15 16]
C =[1 1 ; 2 2]

A*C

A.*B

%Unvectorized vs. Vectorized Implementation
theta = [1;1;2;3];
x = [1;1;2;3];
n = length(x)-1;

prediction = 0.0;
for j = 1:n+1
    prediction = prediction + theta(j)*x(j);
end
disp(prediction)

predictionVect = theta'*x


A = [1 2; 3 4; 5 6]
B = [1 2 3; 4 5 6]

C = magic(10)
D = [1:10]'

C*D
%D'*C
v = zeros(10, 1)
for i = 1:10
    for j = 1:10
        v(i) = v(i) + C(i,j)*D(j);
    end
end
disp(v)


