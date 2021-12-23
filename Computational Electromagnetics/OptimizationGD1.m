%% Advanced Optimization Technique

close all;
clc;
clear all;
%% Settings
options = optimset('GradObj','on', 'MaxIter', 100);
initialTheta = zeros(2,1);
[optTheta, functionVal, exitFlag]...
    =fminunc(@costFunction, initialTheta, options);

disp(optTheta);

%% Cost Function

function [jVal, gradient] = costFunction(theta)

jVal = (theta(1) - 7)^2 + (theta(2) - 6)^2;

%jVal = 3.*(1-theta(1)).^2.*exp(-theta(1).^2 - (theta(2)+1).^2)...
%          - 10.*(theta(1)./5 - theta(1).^3 - theta(2).^5).*exp(-theta(1).^2 - theta(2).^2)...
%            -1/3.*exp(-(theta(1)+1).^2 - theta(2).^2);

gradient = zeros();

gradient(1) = 2*(theta(1) - 7);
gradient(2) = 2*(theta(2) - 6);

%gradient(1) = (2*(theta(1)+1)*exp(-(theta(1)+1)^2-theta(2)^2))/3+6*(1-theta(1))^2*theta(1)*exp(theta(1)^2-(theta(2)+1)^2)-6*(1-theta(1))*exp(theta(1)^2-(theta(2)+1)^2)+20*theta(1)*(-theta(1)^3+theta(1)/5-theta(2)^5)*exp(-theta(1)^2-theta(2)^2)-10*(1/5-3*theta(1)^2)*exp(-theta(1)^2-theta(2)^2);
%gradient(2) = -6*(1-theta(1))^2*(theta(2)+1)*exp(theta(1)^2-(theta(2)+1)^2)+(2*theta(2)*exp(-theta(2)^2-(theta(1)+1)^2))/3+20*theta(2)*(-theta(2)^5-theta(1)^3+theta(1)/5)*exp(-theta(2)^2-theta(1)^2)+50*theta(2)^4*exp(-theta(2)^2-theta(1)^2);

end