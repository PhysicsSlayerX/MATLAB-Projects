%% Reflectance 3D Plot

close all;
clc;
clear all;


S11RotFreq = readtable('ReflectanceRotFreq.csv','PreserveVariableNames',true);
rot_data = S11RotFreq{:,1};
freq_data = S11RotFreq{:,2};
S11_data = S11RotFreq{:,3};
%S11_data = str2double(S11RotFreq{:,3});

S11_matrix = [rot_data freq_data S11_data];
real_S11_matrix = real(S11_matrix);

freq = freq_data(1:50);
rotations = [1:37];

S11rawData = zeros(50, 37);
freqNumber = length(freq_data);


x = 0;
for index = 1:37
    for subIndex = 1+x:50+x
        S11rawData(subIndex-x, index) = S11_data(subIndex);
    end
    x = x + 50;
end

realS11data = real(S11rawData);

fig1 = figure('Color','w');
set(fig1, 'Position', [1 41 1920 963]);
set(fig1, 'Name', 'S11 vs. Rot. vs Freq Grid');
set(fig1, 'NumberTitle', 'off');


subplot(1,2,1);
h = imagesc(rotations, freq, realS11data); hold on;
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('S11 vs. Rotation vs. Frequency');
axis([min(rotations*10) max(rotations*10) min(freq) max(freq)]);
axis tight;
colorbar;
hold off;


subplot(1,2,2);

figs = surf(rotations,freq, realS11data,'FaceAlpha',0.5); hold on;
ylabel('$\textrm{Frequency (GHz)}$','Interpreter','LaTex', 'FontSize', 14);
xlabel('$\textrm{Rotation (deg)}$','Interpreter','Latex', 'FontSize', 14);
zlabel('$\textrm{Re(S11)}$','Interpreter','Latex', 'Rotation', 90, 'FontSize', 14);
title('$\textrm{3D View, Surface Plot of the S11 vs. Rot. vs. Freq}$','Interpreter','Latex');
xlim([min(rotations) max(rotations)]);
ylim([min(freq) max(freq)]);
zlim([min(-0.2) max(1)]);
figs.EdgeColor = 'none';
colorbar;
disp(max(realS11data));
hold off;