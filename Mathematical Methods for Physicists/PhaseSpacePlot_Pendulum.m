%Phase Space Plot for a Pendulum

close all;
clc;
clear all;

%Physical Constants

THETA_0 = pi/3;
THETA_DOT_0 = 0;

thetaList = zeros();
thetadotList = zeros();
% 
for T = [1:1:5000]
    [theta_i, theta_dot] = theta(0.01*T, THETA_0 , THETA_DOT_0);
    thetaList(T) = theta_i;
    thetadotList(T) = theta_dot;
    THETA_0 = theta_i;
    THETA_DOT_0 = theta_dot;
end

%Plot

fig1 = figure('Color','w');
set(fig1, 'Position', [1 41 1920 963]);
set(fig1, 'Name', 'Phase Space Analysis');
set(fig1, 'NumberTitle', 'off');
plot(thetaList, thetadotList,'-b', 'LineWidth', 3); 
hold off;
axis([min(thetaList) max(thetaList) min(thetadotList) max(thetadotList)]);
xlabel('Angle (rad)');
ylabel('Angular Velocity (rad/s)','Rotation',90);
title(['Phase Space']);
h = legend('Phase Plot','Location','NorthEastOutside');
set(h,'LineWidth',2);

function [theta_double_dot] =  get_theta_double_dot(theta, theta_dot)

    g = 9.8;
    L = 2; 
    u = 0.1;
    theta_double_dot = -u*theta_dot - (g/L)*sin(theta);
end

function [theta ,theta_dot] = theta(t, theta, theta_dot)
    delta_t = 0.01;
    for time = [0, delta_t, t]
        
        theta_double_dot = get_theta_double_dot(theta, theta_dot);
        theta = theta + theta_dot*delta_t;
        theta_dot = theta_dot + theta_double_dot*delta_t;
    end
end

