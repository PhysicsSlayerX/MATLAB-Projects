%Velocity vs Time Graph

close all;
clc;
clear all;

%Physical Constants

VEL_0 = 0;


velList = vel(10, VEL_0);
velList2 = 2*velList;
velList3 = 3*velList;


solutionVel = zeros();
for timeRange = [1:1:10]
    solutionVel(timeRange) = 39.2 - 39.2*exp(-timeRange/4);
end

%Plot


fig1 = figure('Color','w');
set(fig1, 'Position', [2 42 958 954]);
set(fig1, 'Name', 'v-t Analysis');
set(fig1, 'NumberTitle', 'off');

Timerange = [0:0.001:10];


plot(Timerange, velList,'-b', 'LineWidth', 3);  hold on;
plot(Timerange, velList2,'-.r', 'LineWidth', 3); hold on; 
plot(Timerange, velList3,'--g', 'LineWidth', 3); hold on;
grid on;
hold off;
axis([min(Timerange) max(Timerange) 0 100]);
xlabel('Time (s)');
ylabel('Velocity','Rotation',90);
title(['Velocity vs. Time Graph']);
h = legend('Integral Curve','Location','NorthEastOutside');
set(h,'LineWidth',2);
text(6,25, '$$ v = \int \left( g - \frac{\gamma}{m}\right) dt$$','Interpreter', 'latex', 'Fontsize', 18,'Color','b')
text(6,55, '$$ v = 2\int \left( g - \frac{\gamma}{m}\right) dt$$','Interpreter', 'latex', 'Fontsize', 18,'Color','r')
text(6,85, '$$ v = 3\int \left( g - \frac{\gamma}{m}\right) dt$$','Interpreter', 'latex', 'Fontsize', 18,'Color','g')

function vel_dot =  get_vel_dot(vel)

    g = 9.8;
    gamma = 1.0;
    m = 4.0;
    vel_dot = g - [(gamma/m)*vel];
end
%Solving the integral curve
function [velList] = vel(t, vel_value)
    delta_t = 0.001;    %Resolution Size
    velList = zeros();
    time = [0: delta_t: t];
    for index = 1:length(time)
        v_dot = get_vel_dot(vel_value);
        vel_value = vel_value + v_dot*delta_t;
        velList(index) = vel_value;
    end
        
end

