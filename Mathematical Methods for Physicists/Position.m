x0 = 2;
v0 = -8;
vf = 0;
t = [0:1:8];
a = 1; %1 m/s^2

x = x0 + v0*t + (1/2)*a*t.^2;
disp(x);