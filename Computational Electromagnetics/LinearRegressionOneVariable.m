%% Linear Regression, Hypothesis, Cost/Objective Function,

close all;
clc;
clear all;


%% Create a Training Set

data = load('ex1data1.txt'); % read comma separated data
X = data(:, 1); y = data(:, 2);

%% Parameters

m = length(X); % number of training examples
X = [ones(m,1),data(:,1)]; % Add a column of ones to x
theta = zeros(2, 1); % initialize fitting parameters
iterations = 1500;
alpha = 0.01;

figure;
set(figure,'Color','k');
plot(X(:,2), y, 'rx', 'MarkerSize', 10); % Plot the data
hold on;
ylabel('Profit in $10,000s'); % Set the y-axis label
xlabel('Population of City in 10,000s'); % Set the x-axis label
hold off;
plot_darkmode;

disp(computeCost(X, y, theta));
disp(computeCost(X, y, [-1;2]));

% Run gradient descent:
% Compute theta
theta = gradientDescent(X, y, theta, alpha, iterations);

% Print theta to screen
% Display gradient descent's result
fprintf('Theta computed from gradient descent:\n%f,\n%f',theta(1),theta(2))

% Plot the linear fit
hold on; % keep previous plot visible
plot(X(:,2), X*theta, '-')
legend('Training data', 'Linear regression')
hold off % don't overlay any more plots on this figure
plot_darkmode;
% Predict values for population sizes of 35,000 and 70,000
predict1 = [1, 3.5] *theta;
fprintf('For population = 35,000, we predict a profit of %f\n', predict1*10000);

predict2 = [1, 7] * theta;
fprintf('For population = 70,000, we predict a profit of %f\n', predict2*10000);

% Visualizing J(theta_0, theta_1):
% Grid over which we will calculate J
theta0_vals = linspace(-10, 10, 100);
theta1_vals = linspace(-1, 4, 100);

% initialize J_vals to a matrix of 0's
J_vals = zeros(length(theta0_vals), length(theta1_vals));

% Fill out J_vals
for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];    
	  J_vals(i,j) = computeCost(X, y, t);
    end
end

%you can use meshgrid function too!

% Because of the way meshgrids work in the surf command, we need to 
% transpose J_vals before calling surf, or else the axes will be flipped
J_vals = J_vals';

% Surface plot
figure;
set(figure,'Color','k');
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');
plot_darkmode;
% Contour plot
figure;
set(figure,'Color','k');
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 3, 20))
xlabel('\theta_0'); ylabel('\theta_1');
hold on;
plot(theta(1), theta(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
plot_darkmode;
function J = computeCost(X, y, theta)
%COMPUTECOST Compute cost for linear regression
%   J = COMPUTECOST(X, y, theta) computes the cost of using theta as the
%   parameter for linear regression to fit the data points in X and y

% Initialize some useful values
m = length(y); % number of training examples
% You need to return the following variables correctly 
J = 0;
% theta0 = theta(1); %For unvectorized
% theta1 = theta(2); %For unvectorized
% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta
%               You should set J to the cost.
%% Unvectorized:
% for index = 1:m
% 
%     Htheta = theta0.*X(index,1) + theta1.*X(index,2);
%     J = J + (Htheta - y(index)).^2;
%    
% end
% 
% J = (1/(2*m))*J;

%% Semi-Vectorized
for index = 1:m
   Htheta = theta.'*X(index,:)'; 
   %{
    Instruction: 
    a.)Solve for Htheta first. Take the transpose of theta 
        and get the dot product with the transpose of the 
        first row of X matrix (note the X matrix has m column 
        and n+1 columns. The n+1 columns are the x0, x1, etc.)
   
   %} 
   J = J + (Htheta - y(index)).^2;
   %{
    Instruction: 
    b.) Update the value of J for all m values. Take note that
        you need to extract the mth value of the vector y and
        subtract it from Htheta - a scalar value now. 
   
   %} 
end

J = (1/(2*m))*J;
    %{
    Instruction: 
    c.) Don't forget to multiply it with the factor outside
        the summation symbol.
   %} 
% =========================================================================

end

function [theta, J_history] = gradientDescent(X, y, theta, alpha, num_iters)
%GRADIENTDESCENT Performs gradient descent to learn theta
%   theta = GRADIENTDESCENT(X, y, theta, alpha, num_iters) updates theta by 
%   taking num_iters gradient steps with learning rate alpha

% Initialize some useful values
m = length(y); % number of training examples
J_history = zeros(num_iters, 1);

diff_times_x = zeros(m,length(theta));
for iter = 1:num_iters

    % ====================== YOUR CODE HERE ======================
    % Instructions: Perform a single gradient step on the parameter vector
    %               theta. 
    %
    % Hint: While debugging, it can be useful to print out the values
    %       of the cost function (computeCost) and gradient here.
    %

    for index = 1:m
        Htheta = theta.'*X(index,:)';
        for elem = 1:length(theta)
            diff_times_x(index, elem) = (Htheta - y(index)).*X(index,elem);
        end
    end
    %{
    Instruction: 
    a.) You need to have a nested for loop, one for looping on
        all values of m and one for looping for the parameters, in
        this case, theta0 and theta1
    b.) The same case of the cost function, calculate Htheta for
        certain values of m. 
    c.) Proceed in calculating the (Htheta-y)x term. But since were
        going to solve for this for all m values and parameters, let's 
        define a matrix diff_times_x with m rows and n+1 columns.
    %} 
    
    for elem = 1:length(theta)
        theta(elem) = theta(elem) - alpha*(1/m)*sum(diff_times_x(:,elem));
    end
    %{
    Instruction: 
    d.) Update the theta using the calculated diff_times_x matrix, but
        in this case, take the sum of the first column of diff_times_x 
        matrix - this is the summation of all the (Htheta - y)x terms over
        m values. 
    e.) Apply the update for all the parameters, indicated by the elem
        indexing in the for loop.
    %} 
    

    % ============================================================

    % Save the cost J in every iteration    
    J_history(iter) = computeCost(X, y, theta);
    disp(J_history(iter));
end
    
end

