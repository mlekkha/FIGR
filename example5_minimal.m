%======== USER-CONFIGURABLE OPTIONS =======
numGenes      = 2;     	% 2-gene toy model
numNuclei     = 20;     % 20 nuclei (initial coords)
numTimepoints = 21;     % 21 timepoints
timestep      = 0.1;    

grn.Tgg     = [1 -.5 ; -.5 1];   % gene 1 and 2 are both autoactivating!!!!!!!!!!
grn.hg      = [-.25    -.25]';   % trajs are supposed to diverge to
grn.Rg      = [2       1]';      % 4 corners
grn.lambdag = [2       1]';
grn.Dg      = [0       0]';

%======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
xntg        = NaN ([numNuclei,numTimepoints,numRegs]);
xntg(:,1,:) = rand (numNuclei, numRegs) .* grn.Rg' ./ grn.lambdag';
xntg(:,:,numGenes+1:end) = 0.;



%========  =======
global opts;
global ODEopts;
clc;
tt      = (0:numTimepoints-1)' * timestep;
numRegs = numGenes;

gmax = numGenes;  % YLL's useful aliases
nmax = numNuclei;
tmax = numTimepoints;
kmax = numNuclei * numTimepoints;

%======== SET OPTIONS (see README.md for description) ========
% NOTE: infer IS NOT BEING CALLED. 
% THE ONLY RELEVANT OPTIONS ARE ODE_XXX
%
% NOTE: slopethresh, exprthresh, splinesmoothing are in pvxOpts_ngo
% NOTE: Rld_tsafety also should be removed.
pvxOpts_ngo = NaN (numNuclei, numGenes, 3);
pvxOpts_ngo(:,:,1) = 1.00;  % p (spline unsmoothing parameters)
pvxOpts_ngo(:,:,2) = 0.01;  % v (velocity thresholds)
pvxOpts_ngo(:,:,3) = 0.20;  % x (expression thresholds)
opts = struct(  'debug', 0, ...
    'Rld_tsafety', 3, ...       % should eventually ged rid
    'spatialsmoothing', 0.5, ...
    'minborder_expr_ratio', 0.01, ...
    'Rld_method', 'slope', ...
    'synthesisfunction', 'synthesis_sigmoid_sqrt', ...
    'ODEAbsTol', 1e-5, ...      % originally 1e-3
    'ODEsolver', 'ode45', ...
    'pvxOpts_ngo', pvxOpts_ngo, ...
    'lambda', 0.5, ...
    'lm', 'FIGRlogReg'); % glmfit|FIGRlogReg|lassoglm
ODEopts = odeset('AbsTol', opts.ODEAbsTol);



%======== COMPUTE TRAJECTORIES FOR THE TOY MODEL ABOVE ========
[xntg] = computeTrajs (opts, grn, xntg, tt);


%======== FIGURE 1: TRAJECTORIES IN  (x_1, x_2) PLANE ========
close all;
figure ('Position', [0 0 600 600]); hold on;
title ('TRAJECTORIES IN  (x_1, x_2) PLANE');
xlabel ('$x_A$', 'interpreter', 'latex');
ylabel ('$x_B$', 'interpreter', 'latex');
axis ([0 1 0 1]);
plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))');

%======== FIGURE 2: TRAJECTORIES x_1(t) ========
figure ('Position', [600 0 600 600]); hold on;
title ('FIGURE 2: TRAJECTORIES x_1(t)');
xlabel ('$t$', 'interpreter', 'latex');
ylabel ('$x_A$', 'interpreter', 'latex');
for n=1:nmax
    plot (tt, xntg(n,:,1));
end

fprintf ("Trajectories should DIVERGE, not CONVERGE as we see!\n");


return;
