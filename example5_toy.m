%======== USER-CONFIGURABLE OPTIONS =======
numGenes      = 2;     	% hardwired for now
numNuclei     = 20;     % number of "nuclei" or "initial conds"
numTimepoints = 21; 
timestep      = 0.1;

%======== CODE STARTS HERE =======
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
%
%
%OVERRRIDING!
grn.Tgg = [1 -.5 ; -.5 1]; 
grn.hg = [-.25 -.25]';
grn.Rg = [2 1]';
grn.lambdag = [2 1]';
%======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
xntg        = NaN ([numNuclei,numTimepoints,numRegs]);
xntg(:,1,:) = rand (numNuclei, numRegs) .* grn.Rg' ./ grn.lambdag';
xntg(:,:,numGenes+1:end) = 0.;
%======== COMPUTE TRAJECTORIES FOR THE TOY MODEL ABOVE ========
[xntg] = computeTrajs (opts, grn, xntg, tt);

%
%
% THERE IS SOMETHING WRONG.  THE ABOVE CODE
% NEVER GIVES TRAJS THAT DIVERGE (EVEN IF GENE 1 is autoactivatin!)
%
% TO DO: PUT IN A T MATRIX AND OTHER PARS
% THAT ARE KNOWN TO DIVERGE TO 4 CORNERS.
%   SEE IF WE STILL SEE CONVERGENCE...
disp (grn.Tgg); disp (grn.hg); disp (grn.Rg); disp (grn.lambdag);





%======== RUN infer
[grnFIGR, diagnostics] = infer (opts, xntg, tt, numGenes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======== VISUALIZE GENE STATES ========
close all;
figure ('Position', [0 0 600 600]); hold on;
set (gcf, 'name', 'Data');
title ('Gene A velocity (green $v_A>0$, red $v_A<0$)', 'interpreter','latex');
xlabel ('$x_A$', 'interpreter', 'latex');
ylabel ('$x_B$', 'interpreter', 'latex');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.4 .4 .4]); %gray
%======== Datapoints (filled circles) ========
gTarget = 1;


xkg = reshape (xntg, [kmax gmax]);
%ykg = reshape (diagnostics.yntg, [kmax gmax]);
%vkg = reshape (diagnostics.vntg, [kmax gmax]);
yk = diagnostics.yntg(:,:,gTarget); yk = yk (:);
vk = diagnostics.vntg(:,:,gTarget); vk = vk (:);

markerSizes  = 100.0 * abs (vk) .^ 0.5;  % bigger markers mean stronger vels
markerColors = NaN (kmax, 3);
for k=1:kmax
    if (yk(k)>0) 
        markerColors(k,:) = [0 .5 0]; % green = ON
    elseif (yk(k)<0) 
        markerColors(k,:) = [.7 0 0]; % red = OFF
    else
        markerColors(k,:) = [.5 .5 .5]; % gray = UNKNOWN
    end
end
scatter (xkg(:,1), xkg(:,2), markerSizes, markerColors, 'filled');


return;




%======== COMPUTE MODEL PREDICTIONS FOR VELOCITIES (fitted v's) ========
% Use the MODEL (grnFIGR) to PREDICT the values of v at the ORIGINAL xkg points.
Tgg = grnFIGR.Tgg;
hg = grnFIGR.hg;
Rg = grnFIGR.Rg;
lambdag = grnFIGR.lambdag;
% Tgg = grn.Tgg;
% hg = grn.hg;
% Rg = grn.Rg;
% lambdag = grn.lambdag;
ykg = NaN (kmax, gmax);   % this is DIFFERENT from earlier!
vkg = NaN (kmax, gmax);   % these 
for g=1:gmax
    for k=1:kmax
        vkg(k,g) = ...
            Rg(g) * heaviside( Tgg(g,:) * xkg(k,:)'  + hg(g) ) ...
            - lambdag(g) * xkg(k,g);
    end
end
ykg = sign (vkg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure ('Position', [600 0 600 600]); hold on;
set (gcf, 'name', 'Model (fit)');
title ('Gene A velocity (green $v_A>0$, red $v_A<0$)', 'interpreter','latex');
xlabel ('$x_A$', 'interpreter', 'latex');
ylabel ('$x_B$', 'interpreter', 'latex');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
%plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.4 .4 .4]); %gray
%======== Datapoints (filled circles) ========
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

%======== RECOMPUTE TRAJECTORIES
%[xntgRECAL] = computeTrajs (opts, grnFIGR, xntg, tt);

%======== RERUN infer ???
%[grnDUMMY, diagnosticsRECAL] = infer (opts, xntgRECAL, tt, numGenes);
