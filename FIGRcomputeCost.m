function [J,grad] = FIGRcomputeCost(Theta,X,y,lambda)

% This function computes the cost function for regularized logistic

% regression. J is the cost function with uregularized part

%           (1 / m) * sum( - y .* log(h) - (1-y) .* log(1 - h) )

% and regularization term

%           (lambda / (2 * m)) * sum(Theta(2:end).^2).

% grad(i) is the partial derivative of the cost function J with respect to

% Theta parameters. Theta parameters correspond to T and h parameters in the 

% gene circtui model. Here we are not regularizing the bias term.

m = size(X,1);

h = sigmoid(X * Theta);

%% Regularization
J = (1 / m) * sum( - y .* log(h) - (1-y) .* log(1 - h) ) + (lambda / (2 * m)) * sum(Theta(2:end).^2);

grad = zeros(size(Theta, 1), 1);

grad(1) = sum(X(:,1).*(h - y))/m;

for i = 2 : size(grad)
    grad(i) = (1 / m) * sum( (h - y)' * X(:, i) ) + (lambda / m) * sum(Theta(2:end).^2);
end

end