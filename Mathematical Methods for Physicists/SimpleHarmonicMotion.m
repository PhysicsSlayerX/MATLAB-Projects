

clc;
close all;
clear all;

POS_X = 10; %Starting position - Amplitude
X_DOT = 0; %Initial Velocity
T = 15;

[x,xd] = get_pos_x(T, POS_X, X_DOT);

%Plot
Timerange = [0:0.01:T];

fig1 = figure('Color','w');
set(fig1, 'Position', [2 42 958 954]);
set(fig1, 'Name', 'Position vs. Time Analysis');
set(fig1, 'NumberTitle', 'off');

subplot(1, 2, 1);
plot(Timerange, x,'-b', 'LineWidth', 3);  hold on;
plot(Timerange, xd,'-.r', 'LineWidth', 3); hold on; 
grid on;
hold off;
axis([min(Timerange) max(Timerange) -11 11]);
xlabel('Time (s)');
ylabel('Position','Rotation',90);
title(['Position vs. Time Graph']);
h = legend('Postion','Velocity','Location','NorthEastOutside');
set(h,'LineWidth',2);

subplot(1, 2, 2);
plot(x, xd,'-g', 'LineWidth', 3);  hold on; 
grid on;
hold off;
axis([-11 11 -11 11]);
xlabel('Position(m)');
ylabel('Velocity','Rotation',90);
title(['Phase Diagram - Simple Harmonic Motion']);
h = legend('Postion','Velocity','Location','NorthEastOutside');
set(h,'LineWidth',2);


function [x_double_dot] = get_x_double_dot(pos_x, x_dot)
    k = 1;   %spring constant of the oscillator
    m = 1;   %mass of the oscillator
    b = 0;   %damping coefficient
    x_double_dot = -b*x_dot - ((k/m)^2)*pos_x;

end

function [pos_x_List, x_dot_List] = get_pos_x(t, pos_x, x_dot)
    delta_t = 0.01;
    timerange = [0:delta_t:t];
    pos_x_List = zeros();
    x_dot_List = zeros();
    for time = 1:length(timerange)
        x_double_dot = get_x_double_dot(pos_x, x_dot);
        pos_x = pos_x + x_dot*delta_t;
        x_dot = x_dot + x_double_dot*delta_t;
        pos_x_List(time) = pos_x;
        x_dot_List(time) = x_dot;
    end
end
