%==============================================
% VARIABLE NAMING CONVENTIONS IN THIS SCRIPT
%
%  k=1:kmax    data index
%  g=1:gmax    component index
%  x           independent variable
%  y           dependent variable
%  b           model parameter (beta)
%
%==============================================

xkg = [1 2 3 4 5 6 7 8 9]';
yk  = [0 0 0 0 1 0 1 1 1]';

clc; close all; tiledlayout ('flow');

%======== USE MATLAB BUILTIN glmfit()
bg = glmfit (xkg, yk, 'binomial', 'link', 'logit');
ykFit = glmval (bg, xkg, 'logit');

fprintf ("Optimized model parameters from glmfit():  bg = \n") ; disp (bg);
nexttile;
title ('glmfit()'); xlabel ('x'); ylabel ('y');
axis ([0 10 -0.1 1.1]); hold on;  
plot (xkg, yk, 'bo');
plot (xkg, ykFit, 'r*');
xkgDense = 0:.01:10; plot (xkgDense, glmval (bg, xkgDense, 'logit'), 'r-');

%======== USE HOMEMADE maximizeLogL()
bg = maximizeLogL (xkg, yk);
ykFit = computeModelPrediction (bg, xkg);

fprintf ("Optimized model parameters from maximizeLogL():  bg = \n") ; disp (bg);
nexttile;
title ('maximizeLogL()'); xlabel ('x'); ylabel ('y');
axis ([0 10 -0.1 1.1]); hold on;  
plot (xkg, yk, 'bo');
plot (xkg, ykFit, 'r*');
xkgDense = 0:.01:10; plot (xkgDense, glmval (bg, xkgDense, 'logit'), 'r-');
return;

%=====================================================
%  COMPUTE MEAN FUNCTION AT EACH DATA POINT
%=====================================================
function yk = computeModelPrediction (bg, xkg)
kmax = size (xkg,1);
xkg = [ones(kmax,1)  xkg];    % prepend column of 1's
yk = 1 ./ (1 + exp(-xkg*bg));
end
%=====================================================
%  COMPUTE LOG LIKELIHOOD OF PARAMETER VECTOR bg
%=====================================================
function logL = computeLogL (bg, xkg, yk)
kmax = size (xkg,1);
xkg = [ones(kmax,1)  xkg];    % prepend column of 1's
yk = 2*yk - 1;                % convert to y=+1/-1 convention
logL = -sum( log (1 + exp( -yk .* (xkg*bg)  )) );
end
%=====================================================
%  MAXIMIZE LOG LIKELIHOOD 
%  (PERFORM LOGISTIC REGRESSION)
%=====================================================
function bg = maximizeLogL (xkg, yk)
gmax = size (xkg,2);
bgInitial = zeros (gmax+1, 1);  % zero initial conditions
bg = fminsearch (@(bg) -computeLogL(bg,xkg,yk), bgInitial);
end



% Alternative method:
% b0 = bg(1);      % split parameter array into
% bg = bg(2:end);  % b0 and [b1,b2,...,bg] parts
% logL = sum( log (1 + exp( -yk .* (b0 + xkg*bg)  )) );
%[kmax gmax] = size (xkg);
