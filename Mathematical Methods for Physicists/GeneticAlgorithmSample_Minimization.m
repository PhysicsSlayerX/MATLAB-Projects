%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Genetic Algorithm Code - Maximization Problem 
%% Consider this to be refactored!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%---------------------------------------------
p = 100; %population size
c = 30; % number of pairs of chromosomes to be crossovered
m = 30; %number of chromosomes to be mutated
tg = 250; % Total number of generations

fig1 = figure('Color','w');
set(fig1, 'Position', [1 41 1920 963]);
set(fig1, 'Name', 'Genetic Algorithm');
set(fig1, 'NumberTitle', 'off');

subplot(2,2,[1,3]);

XX = -3:6/100:3;
YY = XX;
[YYY,XXX] = meshgrid(XX);
ZZZ = 3.*(1-XXX).^2.*exp(-XXX.^2 - (YYY+1).^2)...
            - 10.*(XXX./5 - XXX.^3 - YYY.^5).*exp(-XXX.^2 - YYY.^2)...
            -1/3.*exp(-(XXX+1).^2 - YYY.^2);
figs = surf(XXX,YYY,ZZZ.','FaceAlpha',0.5); hold on;
figs1 = get(figs, 'Parent');
set(figs1,'FontSize',10,'LineWidth',0.5);
%set(figs1,'Xdir','reverse')
xlabel('$y$','Interpreter','LaTex', 'FontSize', 14);
ylabel('$x$','Interpreter','Latex', 'FontSize', 14);
zlabel('$z$','Interpreter','Latex', 'Rotation', 90, 'FontSize', 14);
title('$\textrm{3D View, Surface Plot}$','Interpreter','Latex');
xlim([-3 3]);
ylim([-3 3]);
zlim([-10 10]);
figs.EdgeColor = 'none';
colorbar;
%disp(max(ZZZ));
hold off;

subplot(2,2,2);
imagesc(XX, YY, ZZZ.'); hold on;
h = imagesc(XX,YY, ZZZ.');
h2 = get(h, 'Parent');
set(h2,'FontSize',10,'LineWidth',0.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','right');
title('$\textrm{Top View, Surface Plot}$','Interpreter','Latex');
axis([-1 +1 -1 +1]);
axis equal tight;
colorbar;
hold off;

subplot(2,2,4);

P = population(p);
K = 0;
[x1, y1] = size(P);
P1 = 0;



for i = 1:tg
    Cr = crossover(P,c);
    Mu = mutation(P,m);
    P(p+1:p+2*c,:) = Cr;
    P(p+2*c+1:p+2*c+m,:) = Mu;
    E = evaluation(P);
    [P, S] = selection(P,E,p);
    K(i,1) = sum(10^6-S)/p;
    K(i,2) = 10^6-S(1); %best
    visuals1 = plot(K(:,1), 'b-o', 'linewidth', 2); 
    hold on;
    visuals2 = plot(K(:,2), 'r-.','linewidth', 2); 
    title('$ \textrm{Genetic Algorithm - Minimization Problem} $','Interpreter','latex','FontSize', 18);
    axis([0 tg -10 10]);  
    xlabel('$ \textrm{Number of Generations} $','Interpreter','latex');
    ylabel('$ \textrm{Minimum Value} $','Interpreter','latex');
    h = legend('Average Value','Minimum','Location','NorthEast');
    set(h,'LineWidth',2);
    hold off;
    grid on;
    drawnow;
 
end    

Min_fitness_value = min(K(:,2));
P2 = P(1,:); %Best Chromosomes
%Convert binary to real numbers
A = bi2de(P2(1,1:y1/2));     %convert binary vectors to decimal
x = -3+A*(3-(-3))/(2^(y1/2) - 1);
B = bi2de(P2(1,y1/2 + 1:y1));
y = -3+B*(3-(-3))/(2^(y1/2) - 1);
Optimal_Solution = [x y];
disp(Optimal_Solution);
disp(Min_fitness_value);

function Y = population(n)
% n - population size
% It is noted that the number of bits to represent the variables 
% in binary numbers depends on the required accuracy (the number of digits
% after comma)

% In this example, I want the solution precision with 5 places after the
% decimal point, and with the upper and lower bounds of the variables are 3
% and -3, so, for each variable, we need 20 bits.

% General formula: 2^(m-1) (upper bound  - lower bound)*10^p < 2*m -1
% In this case: p = 5 and m = 20.

% We have variables (x and y), so we need 40 bits in total for binary
% encoding.

Y = round(rand(n,40));
end

function Y_c = crossover(P,n)
% P = Population
% n = number of pairs of chromosomes to be crossovered

[x1,y1] = size(P);
Z = zeros(2*n, y1);
    for i = 1:n
        r1 = randi(x1, 1, 2);
        while r1(1)== r1(2)
            r1 = randi(x1,1,2);
        end
        A1 = P(r1(1),:); %Parent 1
        A2 = P(r1(2),:); %Parent 2
        r2 = 1 + randi(y1 - 1); %random cutting point
        B1 = A1(1,r2:y1); 
        A1(1,r2:y1) = A2(1,r2:40);
        A2(1, r2:40) = B1;
        Z(2*i-1,:) = A1; %Offspring 1
        Z(2*i,:) = A2; %Offspring 2
    end
Y_c = Z;
end

function Y_m = mutation(P,n)
[x1, y1] = size(P);
Z = zeros(n,y1);
    for i = 1:n
        r1 = randi(x1);
        A1 = P(r1,:);   %rand
        % parent
        r2 = randi(y1);
        if A1(1,r2) == 1
            A1(1,r2) = 0; %flick the bit
        else 
            A1(1,r2) = 1;
        end
        Z(i,:) = A1;
    end
Y_m = Z;
end

function eval = evaluation(P)
[x1, y1] = size(P);
H = zeros(1,x1);
    for i = 1:x1
        A = bi2de(P(i,1:y1/2));     %convert binary vectors to decimal
        x = -3+A*(3-(-3))/(2^(y1/2) - 1);
        B = bi2de(P(i,y1/2 + 1:y1));
        y = -3+B*(3-(-3))/(2^(y1/2) - 1);
        H(1,i) = 3*(1-x)^2*exp(-x^2 - (y+1)^2)...
            - 10*(x/5 - x^3 - y^5)*exp(-x^2 - y^2)...
            -1/3*exp(-(x+1)^2 - y^2); %Objective Function
    end
eval = 10^6 - H;
end

function[YY1, YY2] = selection(P,F,p)
%P = population, F = fitness value, p = population size
[x, y] = size(P);
Y1 = zeros(p, y);
F = F + 10; % adding 10 to ensure no chromosome has negative fitness
%elite selection
Fn = zeros();
e = 3;

    for i = 1:e
        [r1, c1] = find(F == max(F));
        Y1(i,:) = P(max(c1),:);
        P(max(c1),:) = [];
        Fn(i) = F(max(c1));
        F(:, max(c1)) = [];
    end
D = F/sum(F); %Determine selection probability
E = cumsum(D); % Determine cumulative probability
N = rand(1); %Generatre a vector containing normalized random numbers
d1 = 1;
d2 = e;
    while d2 <= p - e
        if N <= E(d1)
            Y1(d2+1,:) = P(d1,:);
            Fn(d2+1) = F(d1);
            N = rand(1);
            d2 = d2+1;
            d1 = 1;
        else
            d1 = d1 + 1;
        end
    end
YY1 = Y1;
YY2 = Fn-10; %substract 10 to return the original fitness
end