function theta = FIGRlogReg(xkg, y, lambda)

% JEH 2020-6
% This function implements regularized logistic regression with cost
% function defined in FIGRlogRegCost. fminunc is used to find the minimum
% of the cost function.
        
[numDatapoints, numGenes] = size(xkg);

xkg = [ones(numDatapoints,1) xkg];     % /Prepend column of 1's for the bias term

% Set initial parameters to zero.
% Parameters are h, T1, T2, ... represented by a vector beta

initialBeta = zeros(numGenes+1, 1); % Add extra param for the bias term

% YLL 2020-6-2: 'Display' = 'none' --> suppress details of fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400, 'Display', 'none');

theta = fminunc(@(beta)FIGRlogRegCost(beta, xkg, y, lambda), initialBeta, options);
        
end