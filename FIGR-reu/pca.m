%============================================================
% Principal Component Analysis
% Yen Lee Loh 2020-6-6
%
% WARNING:   Shadows built-in MATLAB "pca" function!
%
% USAGE:     [UggPrime,xkgPrime] = pca (xkg)
% ARGUMENTS: xkg(1:K,1:G)      coordinates x_{kg} of each data point
% RETURNS:   UggPrime(1:G,1:G) principal-axis vectors U_{gg'}
%            xkgPrime(1:K,1:G) coordinates x_{kg'} in principal frame
%
% The algorithm is as follows.  Assume w.l.o.g. that <x_g> = 0 for all g.
% Find the covariance matrix C_{g_1 g_2} = <x_{g_1} x_{g_2}>, 
%   computed over all datapoints k.
% Find the eigenvalue matrix Lambda_{ij} = lambda_i delta{ij}
% and eigenvector matrix U_{ij} = u^(j)_i whose columns are the evecs
% such that the eigenequation C U = U Lambda is satisfied.  
% Now, you can convince yourself that x = U x'
% where x' is a column vector whose components are the new coordinates.
% Therefore x' = U^{-1} x.  In index notation, X_{g'g} = U_{gg'} and
%   x_{g'} = X_{g'g} x_{g}.
%
% For the general problem, the principal directions are the columns of U,
%   U_{gg'},
% and the datapoint coordinates in the principal frame are 
%   x' = U^{-1} (x - x0).
%
%============================================================
function [UggPrime, xkgPrime] = pca (xkg)
fprintf ('\n============== pca() =====================\n');

xgMean = mean(xkg,1)';     % CoM over all data points (average over k)
Cgg = cov (xkg);           % covariance matrix

[evecs,evals] = eig (Cgg); % eigenvalues and eigenvectors
evals = diag (evals);
[~,ind] = sort (evals);    % find indices that sort evals in ascending order
ind = flip (ind);
evals = evals(ind);        % sort evals in descending order of evals
evecs = evecs(:,ind);      % sort evecs in descending order of evals
evecs(:,2) = sign(det(evecs))*evecs(:,2); % hack to make det U=+1
UggPrime = evecs;

xgk = xkg';
xgPRIMEk = inv(UggPrime) * (xgk - xgMean);
xkgPrime = xgPRIMEk';  % convert back to correct form for returnvalue

fprintf ("Center of mass:\n"); disp (xgMean);
fprintf ("Covariance matrix:\n"); disp (Cgg);
fprintf ("Eigenvals in decreasing order:\n"); disp (evals'); 
fprintf ("Eigenvecs in decreasing order of evals (columns below):\n"); disp (evecs);
fprintf ("Principal directions (columns below):\n"); disp (UggPrime);
fprintf ("Inverse transformation:\n"); disp (inv(UggPrime));
%fprintf ("Orthonormality check:\n"); disp (evecs * evecs');
fprintf ("============== end pca() =====================\n");
return
end
