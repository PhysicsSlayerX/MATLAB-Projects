%% Logistic Regression

close all;
clc;
clear all;

% Load Data
% The first two columns contain the exam scores and the third column contains the label.
data = load('ex2data1.txt');
X = data(:, [1, 2]); 
y = data(:, 3);

% Plot the data with + indicating (y = 1) examples and o indicating (y = 0) examples.
plotData(X, y);

%  Setup the data matrix appropriately
[m, n] = size(X);

% Add intercept term to X
X = [ones(m, 1) X];

% Initialize the fitting parameters
initial_theta = zeros(n + 1, 1);

% Compute and display the initial cost and gradient
[cost, grad] = costFunction(initial_theta, X, y);
fprintf('Cost at initial theta (zeros): %f\n', cost);
disp('Gradient at non-zero theta:'); disp(grad);

% Compute and display cost and gradient with non-zero theta
test_theta = [-24; 0.2; 0.2];
[cost, grad] = costFunction(test_theta, X, y);
fprintf('\nCost at non-zero test theta: %f\n', cost);
disp('Gradient at non-zero theta:'); disp(grad);

%  Set options for fminunc
options = optimoptions(@fminunc,'Algorithm','Quasi-Newton','GradObj', 'on', 'MaxIter', 400);

%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[theta, cost] = fminunc(@(t)(costFunction(t, X, y)), initial_theta, options);

% Print theta
fprintf('Cost at theta found by fminunc: %f\n', cost);
disp('theta:');disp(theta);
% Plot Boundary
plotDecisionBoundary(theta, X, y);
% Add some labels 
hold on;
% Labels and Legend
xlabel('Exam 1 score')
ylabel('Exam 2 score')
% Specified in plot order
legend('Admitted', 'Not admitted')
hold off;

%% Evaluation of Logistic Regression

%  Predict probability for a student with score 45 on exam 1  and score 85 on exam 2 
prob = sigmoid([1 45 85] * theta);
fprintf('For a student with scores 45 and 85, we predict an admission probability of %f\n\n', prob);
% Compute accuracy on our training set
p = predict(theta, X);
fprintf('Train Accuracy: %f\n', mean(double(p == y)) * 100);


%% Regularized Logistic Regression

%  The first two columns contains the X values and the third column
%  contains the label (y).
data = load('ex2data2.txt');
XX = data(:, [1, 2]); yy = data(:, 3);

plotData(XX, yy);
% Put some labels 
hold on;
% Labels and Legend
xlabel('Microchip Test 1')
ylabel('Microchip Test 2')
% Specified in plot order
legend('y = 1', 'y = 0')
hold off;

% Add Polynomial Features
% Note that mapFeature also adds a column of ones for us, so the intercept term is handled
XX = mapFeature(XX(:,1), XX(:,2));

% Initialize fitting parameters
initial_theta = zeros(size(XX, 2), 1);

% Set regularization parameter lambda to 1
lambda = 1;

% Compute and display initial cost and gradient for regularized logistic regression
[cost, grad] = costFunctionReg(initial_theta, XX, yy, lambda);
fprintf('Cost at initial theta (zeros): %f\n', cost);
fprintf('Expected cost (approx): 0.693\n');
fprintf('Gradient at initial theta (zeros) - first five values only:\n');
fprintf(' %f \n', grad(1:5));
fprintf('Expected gradients (approx) - first five values only:\n');
fprintf(' 0.0085\n 0.0188\n 0.0001\n 0.0503\n 0.0115\n');



% Compute and display cost and gradient with all-ones theta and lambda = 10
test_theta = ones(size(XX,2),1);
[cost, grad] = costFunctionReg(test_theta, XX, yy, 10);
fprintf('\nCost at test theta (with lambda = 10): %f\n', cost);
fprintf('Expected cost (approx): 3.16\n');
fprintf('Gradient at test theta - first five values only:\n');
fprintf(' %f \n', grad(1:5));
fprintf('Expected gradients (approx) - first five values only:\n');
fprintf(' 0.3460\n 0.1614\n 0.1948\n 0.2269\n 0.0922\n');

%Using fminunc
% Initialize fitting parameters
initial_theta = zeros(size(XX, 2), 1);

lambda = 1;
% Set Options
options = optimoptions(@fminunc,'Algorithm','Quasi-Newton','GradObj', 'on', 'MaxIter', 1000);

% Optimize
[theta, J, exit_flag] = fminunc(@(t)(costFunctionReg(t, XX, yy, lambda)), initial_theta, options);

% Plot Boundary
plotDecisionBoundary(theta, XX, yy);
hold on;
title(sprintf('lambda = %g', lambda))

% Labels and Legend
xlabel('Microchip Test 1')
ylabel('Microchip Test 2')

legend('y = 1', 'y = 0', 'Decision boundary')
hold off;

% Compute accuracy on our training set
p = predict(theta, XX);


fprintf('Train Accuracy: %f\n', mean(double(p == yy)) * 100);


%% Functions

function plotData(X, y)
%PLOTDATA Plots the data points X and y into a new figure 
%   PLOTDATA(x,y) plots the data points with + for the positive examples
%   and o for the negative examples. X is assumed to be a Mx2 matrix.

% Create New Figure
figure; hold on;
set(figure, 'Color', 'k');

% ====================== YOUR CODE HERE ======================
% Instructions: Plot the positive and negative examples on a
%               2D plot, using the option 'k+' for the positive
%               examples and 'ko' for the negative examples.
%
% Plot the data with + indicating (y = 1) examples and o indicating (y = 0) examples.

% Find Indices of Positive and Negative Examples
pos = find(y==1); neg = find(y == 0);
% Plot Examples
plot(X(pos, 1), X(pos, 2), 'k+','LineWidth', 2, 'MarkerSize', 7); hold on
plot(X(neg, 1), X(neg, 2), 'ko', 'MarkerFaceColor', 'y','MarkerSize', 7);
hold off;
% Labels and Legend
xlabel('Exam 1 score')
ylabel('Exam 2 score')

% Specified in plot order
legend('Admitted', 'Not admitted')
plot_darkmode;


% =========================================================================



hold off;

end


function g = sigmoid(z)
%SIGMOID Compute sigmoid function
%   g = SIGMOID(z) computes the sigmoid of z.

% You need to return the following variables correctly 
g = zeros(size(z));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the sigmoid of each value of z (z can be a matrix,
%               vector or scalar).

%Unvectorized
for index = 1:size(z,1)
    for index2 = 1:size(z,2)
    g(index,index2)= 1/(1 + exp(-z(index,index2)));
    end
end

%{
Instructions:
For the sigmoid function to work on any matrices, vectors,
or scalars, we must assume first that z is a matrix.
getting the size of z will make us know its dimension.

We used the for loops for all rows and column elements
to apply the sigmoid function, and store it in "g" with
the same size of the z matrix input. 
%}


% =============================================================

end

function [J, grad] = costFunction(theta, X, y)
%COSTFUNCTION Compute cost and gradient for logistic regression
%   J = COSTFUNCTION(theta, X, y) computes the cost of using theta as the
%   parameter for logistic regression and the gradient of the cost
%   w.r.t. to the parameters.

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta.
%               You should set J to the cost.
%               Compute the partial derivatives and set grad to the partial
%               derivatives of the cost w.r.t. each parameter in theta
%
% Note: grad should have the same dimensions as theta
%

%Semi-vectorized
% for mindex = 1:m
%    J = J + -y(mindex).*log(sigmoid(theta.'*X(mindex,:)'))-...
%        (1 - y(mindex)).*log(1 - sigmoid(theta.'*X(mindex,:)'));
%    % or
%    % J = J + -y(mindex).*log(sigmoid(sum(theta.*X(mindex,:)')))-...
%    %    (1 - y(mindex)).*log(1 - sigmoid(sum(theta.*X(mindex,:)')));
% end
% J = 1/m*J;

%{
Instructions:
For the cost function to work, we must understand that it
will spit a single value for all the parameters involved.
this is the function we want to minimized. 

We check the dimensions so that we will end up with a single
value, not an array of numbers. We take theta.'*X(mindex,:)'
and put it in the sigmoid function and this will spit a
singular value. (In this example, we are iterating over
the m values, so we take the row values for matrix X over
all columns, X(mindex,:) - a 1 x n+1 matrix. For the theta
it was initialized as m x 1 matrix. For the dot product
to work, we need to take their transposes (column -> row 
for theta and row to column for X at mindex) and multiply
them matrix wise. Also, we can do an element wise then
use the sum function afterwards - sum(theta.*X(mindex,:)'

Then the iteration will be done all over m values and
then the accumulated J value will be multiplied to 1/m to 
get the final evaluated cost function.
%}


%Vectorized

h = sigmoid(X*theta);
J = (1/m).*(-y.'*log(h) - (1 - y).'*log(1 - h));

%Semi-vectorized
% for nindex = 1:size(theta)
%     for mindex = 1:m
%        grad(nindex) = grad(nindex) + ((sigmoid(theta.'*X(mindex,:)') - y(mindex))*X(mindex,nindex));
%     end
% end
% 
% grad = (1/m)*grad;

%Vectorized
 grad = (1/m)*X.'*(sigmoid(X*theta) - y);

%{
Instruction:

The same case for coding the cost function, we first
get the gradient by taking derivative of the cost function
with respect to theta_j. We need to take the sum for all m values
in the X matrix, and do it for all the parameters in the
theta vector. Finally, multiply it with 1/m to get the
finalized gradient values. 

Take note, this is an semi-vectorized implementation.

%}






% =============================================================

end


function plotDecisionBoundary(theta, X, y)
%PLOTDECISIONBOUNDARY Plots the data points X and y into a new figure with
%the decision boundary defined by theta
%   PLOTDECISIONBOUNDARY(theta, X,y) plots the data points with + for the 
%   positive examples and o for the negative examples. X is assumed to be 
%   a either 
%   1) Mx3 matrix, where the first column is an all-ones column for the 
%      intercept.
%   2) MxN, N>3 matrix, where the first column is all-ones

% Plot Data
plotData(X(:,2:3), y);
hold on

if size(X, 2) <= 3
    % Only need 2 points to define a line, so choose two endpoints
    plot_x = [min(X(:,2))-2,  max(X(:,2))+2];

    % Calculate the decision boundary line
    plot_y = (-1./theta(3)).*(theta(2).*plot_x + theta(1));

    % Plot, and adjust axes for better viewing
    plot(plot_x, plot_y)
    
    % Legend, specific for the exercise
    legend('Admitted', 'Not admitted', 'Decision Boundary')
    axis([30, 100, 30, 100])
else
    % Here is the grid range
    u = linspace(-1, 1.5, 50);
    v = linspace(-1, 1.5, 50);

    z = zeros(length(u), length(v));
    % Evaluate z = theta*x over the grid
    for i = 1:length(u)
        for j = 1:length(v)
            z(i,j) = mapFeature(u(i), v(j))*theta;
        end
    end
    z = z'; % important to transpose z before calling contour

    % Plot z = 0
    % Notice you need to specify the range [0, 0]
    contour(u, v, z, [0, 0], 'LineWidth', 2)
end
hold off

end


function p = predict(theta, X)
%PREDICT Predict whether the label is 0 or 1 using learned logistic 
%regression parameters theta
%   p = PREDICT(theta, X) computes the predictions for X using a 
%   threshold at 0.5 (i.e., if sigmoid(theta'*x) >= 0.5, predict 1)

m = size(X, 1); % Number of training examples

% You need to return the following variables correctly
p = zeros(m, 1);

% ====================== YOUR CODE HERE ======================
% Instructions: Complete the following code to make predictions using
%               your learned logistic regression parameters. 
%               You should set p to a vector of 0's and 1's
%
prob = sigmoid(X*theta);
for mindex = 1:m
    if prob(mindex) >= 0.5
        p(mindex) = 1;
    else
        p(mindex) = 0;
    end  
end

disp(prob);
disp(p);





% =========================================================================


end

function out = mapFeature(X1, X2)
% MAPFEATURE Feature mapping function to polynomial features
%
%   MAPFEATURE(X1, X2) maps the two input features
%   to quadratic features used in the regularization exercise.
%
%   Returns a new feature array with more features, comprising of 
%   X1, X2, X1.^2, X2.^2, X1*X2, X1*X2.^2, etc..
%
%   Inputs X1, X2 must be the same size
%

degree = 6;
out = ones(size(X1(:,1)));
for i = 1:degree
    for j = 0:i
        out(:, end+1) = (X1.^(i-j)).*(X2.^j);
    end
end

end

function [J, grad] = costFunctionReg(theta, X, y, lambda)
%COSTFUNCTIONREG Compute cost and gradient for logistic regression with regularization
%   J = COSTFUNCTIONREG(theta, X, y, lambda) computes the cost of using
%   theta as the parameter for regularized logistic regression and the
%   gradient of the cost w.r.t. to the parameters. 

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta.
%               You should set J to the cost.
%               Compute the partial derivatives and set grad to the partial
%               derivatives of the cost w.r.t. each parameter in theta

%Semi-vectorized
for mindex = 1:m
   J = J + -y(mindex).*log(sigmoid(theta.'*X(mindex,:)'))-...
       (1 - y(mindex)).*log(1 - sigmoid(theta.'*X(mindex,:)'));
   % or
   % J = J + -y(mindex).*log(sigmoid(sum(theta.*X(mindex,:)')))-...
   %    (1 - y(mindex)).*log(1 - sigmoid(sum(theta.*X(mindex,:)')));
end

J = 1/m*J + (lambda/(2*m))*sum(theta(2:1:end).^2);
%Do not regularized theta(1).
%The regularized term is (lambda/(2*m))*sum(theta(2:1:end).^2)

%First Element of the Gradient
for mindex = 1:m
       grad(1) = grad(1) + ((sigmoid(theta.'*X(mindex,:)') - y(mindex))*X(mindex,1));
end
%Semi-vectorized
for nindex = 2:size(theta)
    for mindex = 1:m
       grad(nindex) = grad(nindex) + ((sigmoid(theta.'*X(mindex,:)') - y(mindex))*X(mindex,nindex));
    end
end

grad = (1/m)*grad;
%Retrieve the first element of the grad since it does not
%belong to the regularization.
grad1 = grad(1);
%Regularized theta(2) to the last element (this is
%theta_1 to theta_n.
grad2end = grad(2:1:end) + (lambda/m).*theta(2:1:end);
%Then combine the grad1 and grad2end into a single 
%column vector.
grad = [grad1; grad2end];



% =============================================================

end

