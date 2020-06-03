function theta = FIGRlogReg(X, y, lambda)

% This function implements regularized logistic regression with cost
% function defined in FIGRcomputeCost. fminunc is used to find the minimum
% of the cost function.
        
        [numDatapoints, numGenes] = size(X);
        
        X = [ones(numDatapoints,1) X];     % /Prepend column of 1's
        
        % Set initial parameters to zero.
        % Parameters are h, T1, T2, ...
        
        initialTheta = zeros(numGenes+1, 1); 
              
        % YLL 2020-6-2: 'Display' = 'none' --> suppress details of fminunc
        options = optimset('GradObj', 'on', 'MaxIter', 400, 'Display', 'none');
        
        theta = fminunc(@(Theta)FIGRcomputeCost(Theta, X, y, lambda), initialTheta, options);
        
end