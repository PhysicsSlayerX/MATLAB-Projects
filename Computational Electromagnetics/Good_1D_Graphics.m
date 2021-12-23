% CALCULATE DATA
phi = linspace(0,2*pi,1000);
f1 = cos(phi);
f2 = sin(phi);
% OPEN FIGURE WINDOW
figure('Color','w');
% PLOT DATA
h = plot(phi,f1,'-b','LineWidth',3);
hold on;
plot(phi,f2,'--r','LineWidth',3);
hold off;
% SET AXIS LIMITS
xlim([0 2*pi]);
ylim([-1.1 1.1]);
% MAKE LINES THICK AND FONTS BIGGER
h2 = get(h,'Parent');
set(h2,'LineWidth',3,'FontSize',18);
% SET TICK MARKS
L = {'-1.0' '-0.5' '0' '0.5' '1.0'};
set(h2,'XTick',[0:6],'YTick',[-1:0.5:1],'YTickLabel',L);
% LABEL AXES
xlabel('$ \textrm{Angle, } \theta $','Interpreter','latex');
ylabel('$ f(x) $','Interpreter','latex');
% ADD LEGEND
h = legend('cos(\theta)','sin(\theta)','Location','NorthOutside');
set(h,'LineWidth',2);

