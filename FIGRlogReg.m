function theta = FIGRlogReg(X, y, numNuclei, lambda)

% This function implements regularized logistic regression with cost

% function defined in FIGRcomputeCost. fminunc is used to find the minimum

% of the cost function.
        
        [m, n] = size(X);
        
        X = [ones(m,1) X]; %% Add bias term.
        
        [m, n] = size(X); 
        
        initialTheta = zeros(n, 1); %% Set initial parameters theat to zero.
              
        options = optimset('GradObj', 'on', 'MaxIter', 400);
        theta = fminunc(@(t)FIGRcomputeCost(t, X, y, lambda), initialTheta, options);
        
end