x = linspace(-1,1,100);
y = x.^2;
% Plot Function
fig = figure('Color','w');
h = plot(x,y,'-b','LineWidth',2);
% Set Graphics View
h2 = get(h,'Parent');
set(h2,'FontSize',14,'LineWidth',2);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('FINAL PLOT');
% Set Tick Markings
xm = [-1:0.5:+1];
xt = {};
for m = 1 : length(xm)
xt{m} = num2str(xm(m),'%3.2f');
end
set(h2,'XTick',xm,'XTickLabel',xt);
ym = [0.1:0.1:+1];
yt = {};
for m = 1 : length(ym)
yt{m} = num2str(ym(m),'%2.1f');
end
set(h2,'YTick',ym,'YTickLabel',yt);
% Label Minimum
text(-0.75,0.6,'Cool Curve','Color','b','HorizontalAlignment','left');
text(0,0.03,'min','Color','b','HorizontalAlignment','center');