
close all;
clc;
clear all;

%resonant frequency
meters = 1;
millimeters = meters/1000;
second = 1;
hertz = 1/second;
gigahertz = 1e9*hertz;
c0 = 299792458*meters/second;
pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
e_0 = 8.854e-12;
u_0 = (4*pi)*10^-7;

periodicity = 30*millimeters;
h = 1*millimeters;
w = 3*millimeters;
g = [0.1:0.01:1]*millimeters;
R = 10*millimeters;

ratio_gR = g/R;

Cgap = e_0*(h*w)./g;
Cfrin = e_0*(h+w+g);
Csurf = 2.*e_0.*((h+w)./pi).*log10(4*R./g);
Ctotal = Cgap + Cfrin + Csurf;
L = u_0*(R + (w/2))*(log10((8/(h+w))*(R+(w/2))) - 1/2);
f0 = 1./((2*pi)*sqrt(L*Ctotal));

disp(f0*gigahertz);
%disp(ratio_gR);

fig = figure('Color','w');

set(fig, 'Position', [1 41 1920 963]);
set(fig, 'Name', 'SRR Spectral Response');
set(fig, 'NumberTitle', 'off');
x_axis = ratio_gR;
plot(x_axis, f0,'-b','LineWidth', 3); hold on;

axis([min(x_axis) max(x_axis) min(f0) max(f0)]);
xlabel('g/R');
ylabel('f0[GHz]','Rotation',90);
title('Spectral Response');
h = legend('f0','Location','NorthEastOutside');
set(h,'LineWidth',2);

fig2 = figure('Color','w');

set(fig2, 'Position', [1 41 1920 963]);
set(fig2, 'Name', 'Capacitance');
set(fig2, 'NumberTitle', 'off');

plot(x_axis, Cgap,'-b','LineWidth', 3); hold on;
plot(x_axis, Cfrin,'-g','LineWidth', 3); hold on;
plot(x_axis, Csurf,'-r','LineWidth', 3); hold on;
plot(x_axis, Ctotal,'-k','LineWidth', 3); hold on;

axis([min(x_axis) max(x_axis) min(Ctotal) max(Ctotal)]);
xlabel('g/R');
ylabel('Capacitance','Rotation',90);
title('Capacitance vs. g/R');
h = legend('Cgap','Cfrin','Csurf','Ctotal','Location','NorthEastOutside');
set(h,'LineWidth',2);
