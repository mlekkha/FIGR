%======== USER-CONFIGURABLE OPTIONS =======

numGenes      = 2;     	% hardwired for now
numNuclei     = 20;     % number of "nuclei" or "initial conds"
numTimepoints = 21; 
timestep      = 0.1;

%======== CODE STARTS HERE =======
clc;
tt      = (0:numTimepoints-1)' * timestep;
numRegs = numGenes;

gmax = numGenes;  % YLL's useful aliases
nmax = numNuclei;
tmax = numTimepoints;
kmax = numNuclei * numTimepoints;

%======== SET OPTIONS (see README.md for description) ========
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
    'ODEAbsTol', 1e-3, ...
    'ODEsolver', 'ode45', ...
    'pvxOpts_ngo', pvxOpts_ngo, ...
    'lambda', 0.5, ...
    'lm', 'FIGRlogReg'); % glmfit|FIGRlogReg|lassoglm
ODEopts = odeset('AbsTol', opts.ODEAbsTol);



%======== GENERATE RANDOM REGULATORY PARAMETERS T,h
%	Choose a point r0 within the unit hypercube using rand().
%	Choose a unit random vector T in N=(numGenes+numExternals) dimensions,
%		by generating a vector of N normal deviates and normalizing it.
%	Consider the hyperplane T.(r - r0) = 0, that is, T.r + h = 0
%		where h = -T.r0.
grn.Tgg = NaN (numGenes, numRegs);
grn.hg = NaN (numGenes, 1);
for g=1:numGenes
    r0vec    = rand ([numRegs 1]);
    Tvec     = normrnd (0, 1, [numRegs 1]);
    Tvec     = Tvec ./ norm(Tvec);
    h        = -Tvec' * r0vec;
    grn.Tgg(g,:) = Tvec';
    grn.hg(g)    = h;
end   % The above could be modified so that T is auto-activating or auto-repressing
%======== GENERATE RANDOM KINETIC PARAMETERS R,lambda,D
grn.Rg      = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
grn.lambdag = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
grn.Dg      = zeros (numGenes, 1);
%======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
xntg        = NaN ([numNuclei,numTimepoints,numRegs]);
xntg(:,1,:) = rand (numNuclei, numRegs) .* grn.Rg' ./ grn.lambdag';
xntg(:,:,numGenes+1:end) = 0.;
%======== COMPUTE TRAJECTORIES FOR THE TOY MODEL ABOVE ========
[xntg] = computeTrajs (opts, grn, xntg, tt);




%======== RUN infer
[grnFIGR, diagnostics] = infer (opts, xntg, tt, numGenes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======== VISUALIZE GENE STATES ========
close all;
figure ('Position', [0 0 600 600]); hold on;
set (gcf, 'name', 'Actual trajectories');
title('Velocity of A (green $v>0$, red $v<0$)', 'interpreter','latex');
xlabel ('$x_A$', 'interpreter', 'latex');
ylabel ('$x_B$', 'interpreter', 'latex');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
%plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.4 .4 .4]); %gray
%======== Datapoints (filled circles) ========
xkg = reshape (xntg, [kmax gmax]);
ykg = reshape (diagnostics.yntg, [kmax gmax]);
vkg = reshape (diagnostics.vntg, [kmax gmax]);
gTarget = 1;
markerSizes  = 100.0 * abs (vkg(:,gTarget)) .^ 0.5;  % bigger markers mean stronger vels
markerColors = NaN (kmax, 3);
for k=1:kmax
    if (ykg(k,gTarget)>0) 
        markerColors(k,:) = [0 .5 0]; % green = ON
    elseif (ykg(k,gTarget)<0) 
        markerColors(k,:) = [.7 0 0]; % red = OFF
    else
        markerColors(k,:) = [.5 .5 .5]; % gray = UNKNOWN
    end
end
scatter (xkg(:,1), xkg(:,2), markerSizes, markerColors, 'filled');





return;


% BELOW: DEAD CODE
%
%
% The following is actually NOT what I want to do.
% Rather, I want to take the ORIGINAL points x_kg
% and use the MODEL (grnFIGR) to PREDICT the values of v at these points.
% Needs some coding....
%

%======== RECOMPUTE TRAJECTORIES
[xntgRECAL] = computeTrajs (opts, grnFIGR, xntg, tt);

%======== RERUN infer FOR THE PURPOSES OF GENERATING yntg AND vntg
[grnDUMMY, diagnosticsRECAL] = infer (opts, xntgRECAL, tt, numGenes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======== VISUALIZE GENE STATES ========
figure ('Position', [600 0 600 600]); hold on;
set (gcf, 'name', 'Recomputed trajectories');
title('Velocity of A (green $v>0$, red $v<0$)', 'interpreter','latex');
xlabel ('$x_A$', 'interpreter', 'latex');
ylabel ('$x_B$', 'interpreter', 'latex');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
%plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.4 .4 .4]); %gray
%======== Datapoints (filled circles) ========
xkg = reshape (xntgRECAL, [kmax gmax]);
ykg = reshape (diagnosticsRECAL.yntg, [kmax gmax]);
vkg = reshape (diagnosticsRECAL.vntg, [kmax gmax]);
gTarget = 1;
markerSizes  = 100.0 * abs (vkg(:,gTarget)) .^ 0.5;  % bigger markers mean stronger vels
markerColors = NaN (kmax, 3);
for k=1:kmax
    if (ykg(k,gTarget)>0) 
        markerColors(k,:) = [0 .5 0]; % green = ON
    elseif (ykg(k,gTarget)<0) 
        markerColors(k,:) = [.7 0 0]; % red = OFF
    else
        markerColors(k,:) = [.5 .5 .5]; % gray = UNKNOWN
    end
end
scatter (xkg(:,1), xkg(:,2), markerSizes, markerColors, 'filled');

return;
