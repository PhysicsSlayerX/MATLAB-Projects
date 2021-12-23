clc
clear all;
close all;

syms H p M V x h_bar i;

a = 1;
b = 1;

H = (p^2/(2.*M)) + V;

x = b*[0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0];
H = a*[1 0 0;0 3 0; 0 0 5];

p = x*H - H*x;
disp(p)
