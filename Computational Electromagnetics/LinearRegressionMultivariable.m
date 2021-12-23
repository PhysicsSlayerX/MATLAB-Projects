%% Linear Regression with Multiple Features, Hypothesis, Cost/Objective Function,

close all;
clc;
clear all;


%% Parameters

% Run gradient descent
% Choose some alpha value
alpha = 0.1;
num_iters = 500;



%% Load the Data

data = load('ex1data2.txt');
X = data(:, 1:2);
y = data(:, 3);
m = length(y);

% Print out some data points
% First 10 examples from the dataset
fprintf(' x = [%.0f %.0f], y = %.0f \n', [X(1:10,:) y(1:10,:)]');

% Scale features and set them to zero mean
[X, mu, sigma] = featureNormalize(X);

% Add intercept term to X
X = [ones(m, 1) X];

% Init Theta and Run Gradient Descent 
theta = zeros(3, 1);
[theta, ~] = gradientDescentMulti(X, y, theta, alpha, num_iters);

% Display gradient descent's result
fprintf('Theta computed from gradient descent:\n%f\n%f\n%f',theta(1),theta(2),theta(3));

% Estimate the price of a 1650 sq-ft, 3 br house
% ====================== YOUR CODE HERE ======================
areaNorm = (1650 - mu(1))./sigma(1);
brNorm = (3 - mu(2))./sigma(2);
price = [theta(1) + theta(2)*areaNorm + theta(3)*brNorm]; % Enter your price formula here

% ============================================================

fprintf('\nPredicted price of a 1650 sq-ft, 3 br house (using gradient descent):\n $%f', price);

%%Convergence Comparison:

num_iters = 75;
ALPHA = [1.0 0.3 0.1 0.03 0.01];
markings = ["-b", "-.r", "--g", ":k", "*c" ];

% Init Theta and Run Gradient Descent 
theta = zeros(3, 1);
J_history = zeros(num_iters, length(ALPHA));
for index = 1:length(ALPHA)
    [~, J_history(:,index)] = gradientDescentMulti(X, y, theta, ALPHA(index), num_iters);
end
% Plot the convergence graph
for index = 1:length(ALPHA)
    plot(1:num_iters, J_history(:,index), markings(index), 'LineWidth', 2);
    hold on
    xlabel('Number of iterations');
    ylabel('Cost J');
end
hold off;

%%Normal Equation

% Solve with normal equations:
% Load Data
data = csvread('ex1data2.txt');
X = data(:, 1:2);
y = data(:, 3);
m = length(y);

% Add intercept term to X
X = [ones(m, 1) X];

% Calculate the parameters from the normal equation
theta = normalEqn(X, y);

% Display normal equation's result
fprintf('\nTheta computed from the normal equations:\n%f\n%f\n%f', theta(1),theta(2),theta(3));

% Estimate the price of a 1650 sq-ft, 3 br house. 
% ====================== YOUR CODE HERE ======================
area = 1650;
br = 3;
price = [theta(1) + theta(2)*area + theta(3)*br];

% ============================================================

fprintf('\nPredicted price of a 1650 sq-ft, 3 br house (using normal equations):\n $%f', price);     



%% Necessary Functions

function [X_norm, mu, sigma] = featureNormalize(X)
%FEATURENORMALIZE Normalizes the features in X 
%   FEATURENORMALIZE(X) returns a normalized version of X where
%   the mean value of each feature is 0 and the standard deviation
%   is 1. This is often a good preprocessing step to do when
%   working with learning algorithms.

% You need to set these values correctly

%{
Instructions:

1. Set the X dataset as X_norm. X_norm will be the dataset we're
    going to use for the LR multiple features.
2. We need to get the mean for all the features. The columns of
    X dataset will be the features with m number of elements. 
    Pre-allocate the mu variable with mu = zeros(1, size(X,2)), 
    creating a matrix with 1 row and number of columns the same 
    as the number of X dataset features.
3. Just like for the mean, we need to get the standard deviation
    for each feature. Pre-allocate the sigma with
    sigma = zeros(1, size(X, 2)).
%}
X_norm = X; 
mu = zeros(1, size(X, 2)); %size(X,2) gets the column size of X.
sigma = zeros(1, size(X, 2));

% ====================== YOUR CODE HERE ======================
% Instructions: First, for each feature dimension, compute the mean
%               of the feature and subtract it from the dataset,
%               storing the mean value in mu. Next, compute the 
%               standard deviation of each feature and divide
%               each feature by it's standard deviation, storing
%               the standard deviation in sigma. 
%
%               Note that X is a matrix where each column is a 
%               feature and each row is an example. You need 
%               to perform the normalization separately for 
%               each feature. 
%
% Hint: You might find the 'mean' and 'std' functions useful.
%       

for feature = 1:size(X, 2)
    mu(feature) = mean(X(:,feature));
    sigma(feature) = std(X(:,feature));
    X_norm(:,feature) = (X_norm(:,feature) - mu(feature))./sigma(feature);
end
disp(mu);
%{
Instructions:

1. We need to iterate the normalization for each feature. We use 
    for loop to all of columns to apply the normalization
    to all of the features.
2. Solve for the mean and standard deviation per feature.
3. Normalize all the rows per feature by doing 
    X_i = (X_i - mean_i)./sigma_i, where X_i now is the 
    normalized feature.
%}



% ============================================================

end

function J = computeCostMulti(X, y, theta)
%COMPUTECOSTMULTI Compute cost for linear regression with multiple variables
%   J = COMPUTECOSTMULTI(X, y, theta) computes the cost of using theta as the
%   parameter for linear regression to fit the data points in X and y

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta
%               You should set J to the cost.

%% Semi-Vectorized
% for index = 1:m
%    Htheta = theta.'*X(index,:)'; 
%    %{
%     Instruction: 
%     a.)Solve for Htheta first. Take the transpose of theta 
%         and get the dot product with the transpose of the 
%         first row of X matrix (note the X matrix has m column 
%         and n+1 columns. The n+1 columns are the x0, x1, etc.)
%    
%    %} 
%    J = J + (Htheta - y(index)).^2;
%    %{
%     Instruction: 
%     b.) Update the value of J for all m values. Take note that
%         you need to extract the mth value of the vector y and
%         subtract it from Htheta - a scalar value now. 
%    
%    %} 
% end
% 
% J = (1/(2*m))*J;
    %{
    Instruction: 
    c.) Don't forget to multiply it with the factor outside
        the summation symbol.
   %} 

%% Vectorized

J = (1/(2*m))*((X*theta - y).'*(X*theta - y));

% =========================================================================

end

function [theta, J_history] = gradientDescentMulti(X, y, theta, alpha, num_iters)
%GRADIENTDESCENTMULTI Performs gradient descent to learn theta
%   theta = GRADIENTDESCENTMULTI(x, y, theta, alpha, num_iters) updates theta by
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
    %       of the cost function (computeCostMulti) and gradient here.
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
    J_history(iter) = computeCostMulti(X, y, theta);
    %disp(J_history(iter));
end

end

function [theta] = normalEqn(X, y)
%NORMALEQN Computes the closed-form solution to linear regression 
%   NORMALEQN(X,y) computes the closed-form solution to linear 
%   regression using the normal equations.

theta = zeros(size(X, 2), 1);

% ====================== YOUR CODE HERE ======================
% Instructions: Complete the code to compute the closed form solution
%               to linear regression and put the result in theta.
%

theta = pinv(X'*X)*X'*y;


% ---------------------- Sample Solution ----------------------




% -------------------------------------------------------------


% ============================================================

end

