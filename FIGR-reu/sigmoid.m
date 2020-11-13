function g = sigmoid(z)

% JEH 06-2020
% Logistic function with characteristic sigmoid curve.

    g = 1 ./ (1 + exp(-z));
    
end

