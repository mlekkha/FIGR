%==============================================
% VARIABLE NAMING CONVENTIONS IN THIS SCRIPT
%
%  k=1:kmax    data index
%  g=1:gmax    component index
%  x           independent variable
%  y           dependent variable
%  b           model parameter (beta)
%
% NOTES: IN MATLAB ROUTINES,
%    'binomial' automatically implies ('link', 'logit').
%    We're not using the ('Alpha', alpha) option.
% NOTE:  WE DON'T REGULARIZE THE BIAS TERM
%==============================================

xkg = [1 2 3 4 5 6 7 8 9]';
yk  = [0 0 0 0 1 0 1 1 1]';
[kmax gmax] = size (xkg);

clc; close all; figure ('position', [0 0 600 600]); tiledlayout ('flow');

%======== MATLAB'S glmfit()
bg = glmfit (xkg, yk, 'binomial'); 
ykFit = glmval (bg, xkg, 'logit');

fprintf ("Optimized model parameters from glmfit(lambda=0):  bg = \n") ; disp (bg);

%======== YLL'S maximizeLogL() WITH LAMBDA=0
bg = maximizeLogL (xkg, yk, 0.0);                 %!!!!!!!!!!
ykFit = computeModelPrediction (bg, xkg);

nexttile; axis ([0 10 -0.1 1.1]); hold on;
disp ("maximizeLogL(xkg,yk,lambda=0)") ; disp (bg);
title ("maximizeLogL(xkg,yk,lambda=0)"); xlabel ('x'); ylabel ('y');
plot (xkg, yk, 'bo', xkg, ykFit, 'r*');
xkgDense = 0:.01:10; plot (xkgDense, glmval (bg, xkgDense, 'logit'), 'r-');

%======== YLL'S maximizeLogL() WITH LAMBDA=0.1
bg = maximizeLogL (xkg, yk, 0.1);                %!!!!!!!!!!
ykFit = computeModelPrediction (bg, xkg);

nexttile; axis ([0 10 -0.1 1.1]); hold on;
disp ("maximizeLogL(xkg,yk,lambda=0.1)") ; disp (bg);
title ("maximizeLogL(xkg,yk,lambda=0.1)"); xlabel ('x'); ylabel ('y');
plot (xkg, yk, 'bo', xkg, ykFit, 'r*');
xkgDense = 0:.01:10; plot (xkgDense, glmval (bg, xkgDense, 'logit'), 'r-');

%======== JOANNA'S FIGRlogReg() WITH LAMBDA=0.1
bg = FIGRlogReg (xkg, yk, 0.1);                %!!!!!!!!!!
ykFit = computeModelPrediction (bg, xkg);

nexttile; axis ([0 10 -0.1 1.1]); hold on;
disp ("FIGRlogReg(xkg,yk,lambda=0.1)") ; disp (bg);
title ("FIGRlogReg(xkg,yk,lambda=0.1)"); xlabel ('x'); ylabel ('y');
plot (xkg, yk, 'bo', xkg, ykFit, 'r*');
xkgDense = 0:.01:10; plot (xkgDense, glmval (bg, xkgDense, 'logit'), 'r-');

%======== MATLAB'S lassoglm(): LOGISTIC REGRESSION WITH LASSO REGULARIZATION
[bg1 FitInfo] = lassoglm (xkg, yk, 'binomial', 'Lambda', 0.1/60, 'Alpha', 1e-16);        %!!!!!!!!!!
bg0 = FitInfo.Intercept;
bg = [bg0 bg1]';

disp ("lassoglm(xkg,yk,lambda=0.1)") ; disp (bg);

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
function logL = computeLogL (bg, xkg, yk, lambda)
kmax = size (xkg,1);
xkg = [ones(kmax,1)  xkg];    % prepend column of 1's
yk = 2*yk - 1;                % convert to y=+1/-1 convention
logL = ...
    - sum( log (1 + exp( -yk .* (xkg*bg)  )) )  ...
    - 0.5 * lambda * sum(bg(2:end) .^ 2);
%bg' * bg;
end
%=====================================================
%  MAXIMIZE LOG LIKELIHOOD
%  (PERFORM LOGISTIC REGRESSION)
%=====================================================
function bg = maximizeLogL (xkg, yk, lambda)
gmax = size (xkg,2);
bgInitial = zeros (gmax+1, 1);  % zero initial conditions
bg = fminsearch (@(bg) -computeLogL(bg,xkg,yk,lambda), bgInitial);
end



% Alternative method:
% b0 = bg(1);      % split parameter array into
% bg = bg(2:end);  % b0 and [b1,b2,...,bg] parts
% logL = sum( log (1 + exp( -yk .* (b0 + xkg*bg)  )) );
%[kmax gmax] = size (xkg);
