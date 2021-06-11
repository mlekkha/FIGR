function [J,grad] = FIGRlogRegCost(beta,xkg,y,lambda)

% JEH 2020-6
% This function computes the cost function for regularized logistic
% regression. J is the cost function with unregularized part

%           (1 / numDatapoints) * sum( - y .* log(h) - (1-y) .* log(1 - h) )

% and regularization term

%           (lambda / (2 * numDatapoints)) * sum(beta(2:end).^2).

% grad(i) is the partial derivative of the cost function J with respect to
% beta parameters. Beta parameters correspond to T and h parameters in the 
% gene circuit model. Here the bias term is left out from the regularization.
% 
%
% Manu 07/15/2020: Please note that the objective function below is related
% to the one used by MATLAB's lassoglm() by a factor of 2 and the fact that
% the regularization term is not scaled by the number of datapoints:
%
%           (2 / numDatapoints) * sum( - y .* log(h) - (1-y) .* log(1 - h) )
%
%           + lambda  * sum(beta(2:end).^2).
%
%

numDatapoints = size(xkg,1);

h = sigmoid(xkg * beta);

% Regularization with squared lambdas
J = (1 / numDatapoints) * sum( - y .* log(h) - (1-y) .* log(1 - h) ) + (lambda / (2 * numDatapoints)) * sum(beta(2:end).^2);

% Allocation for partial derivatives
grad = zeros(size(beta, 1), 1);

% Derivative for the bias term
grad(1) = sum(xkg(:,1).*(h - y))/numDatapoints;

% Derivatives for beta_1 to beta_numParams
for i = 2 : size(grad)
    grad(i) = (1 / numDatapoints) * sum( (h - y)' * xkg(:, i) ) + (lambda / numDatapoints) * beta(i);
end

end
