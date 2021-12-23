%% Vector Operations 
close all;
clc;
clear all;

A = [1 1 1];
B = [0 1 1];
C = [0 0 1];
D = [0 -1 -1];
E = [-1 -1 -1];
seed = 1;
rng(seed);

%Prelim Problems
Answer = dot(A, B)*dot(A, cross(B, C))*cross(D,E);
Answer2 = (2.*dot(A, E) + dot(C, E)./2)*cross(cross(D, A), cross(2.*B, C)); 
%disp(Answer2)

Amag = sqrt(sum(A.^2));
Bmag = sqrt(sum(B.^2));
AdotB = acos(dot(A, B)./(Amag.*Bmag));
%disp(AdotB*(180/pi));

Resultant = A + B + C + D + E;
magResultant = sqrt(sum(Resultant.^2));
%disp(magResultant);

F1 = [-20*sin(15) 20*cos(15)];
F2 = [-30*cos(43) -30*sin(43)];

RF = F1 + F2;
EF = -RF;
%disp(EF)

Cvec = [2 -2 3];
d_hat = [1/sqrt(2) 1/sqrt(2) 0];

Cvec_dhat = dot(Cvec, d_hat)*(d_hat) + cross(cross(d_hat,Cvec), d_hat);
%disp(Cvec_dhat);
%R = A + B;

%magR = sqrt(R(1)^2 + R(2)^2);
%thetaR = atan(R(2)/R(1))*(180/pi);


%% Kinematics
x0 = 2;
v0 = -8;
a = 1;

t = linspace(0,8,9);

x = zeros(1,length(t));
v = zeros(1,length(t));

for index = 1:length(t)
    x(index) = x0 + v0.*t(index) + (1./2)*a*t(index).^2;
    v(index) = v0 + a.*t(index); 
end
disp(x);
disp(v);
sum(t)

randMatrix = rand(5,3);
mat = rand(4,6);

tot = [randMatrix(:);mat(:)];
disp(tot);

retrieveMat = reshape(tot(16:39), 4,6);
disp(retrieveMat);

theta = 1;
epsilon = 0.01;

diff = ((2*(theta + epsilon)^3 + 2) - (2*(theta - epsilon)^3 + 2))/2*epsilon;
disp(diff);

