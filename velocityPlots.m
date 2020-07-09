%======================================
% MATLAB SCRIPT
% - Make "velocity plot"
%========================================

%======== USER-EDITABLE PARAMETERS ========
gmax     = 2;      % 2-gene toy model
fmax     = 2;      % EXTERNAL REGULATORS NOT SUPPORTED!
nmax     = 50;     % number of nuclei (initial coords)
tmax     = 11;     % number of timepoints
timestep = 0.1;

grn.Tgg     = [1 -.5 ; -.5 1];   % autoactivating and mutually repressing
grn.hg      = [-.25  ; -.25];    % leads to tetrastable behavior
grn.Rg      = [2     ; 1];
grn.lambdag = [2     ; 1];

% %======== GENERATE RANDOM PARAMETERS T,h
% %	Choose a point r0 within the unit hypercube using rand().
% %	Choose a unit random vector T in N=(numGenes+numExternals) dimensions,
% %		by generating a vector of N normal deviates and normalizing it.
% %	Consider the hyperplane T.(r - r0) = 0, that is, T.r + h = 0
% %		where h = -T.r0.
% grn.Tgg = NaN (numGenes, numRegs);
% grn.hg = NaN (numGenes, 1);
% for g=1:numGenes
%     r0vec    = rand ([numRegs 1]);
%     Tvec     = normrnd (0, 1, [numRegs 1]);
%     Tvec     = Tvec ./ norm(Tvec);
%     h        = -Tvec' * r0vec;
%     grn.Tgg(g,:) = Tvec';
%     grn.hg(g)    = h;
% end   % The above could be modified so that T is auto-activating or auto-repressing
% %======== GENERATE RANDOM PARAMETERS R,lambda
% grn.Rg      = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
% grn.lambdag = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
% grn.Dg      = zeros (numGenes, 1);
% %======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
% xntg        = NaN ([numNuclei,numTimepoints,numRegs]);
% xntg(:,1,:) = rand (numNuclei, numRegs) .* grn.Rg' ./ grn.lambdag';
% xntg(:,:,numGenes+1:end) = 0.;



%======== DERIVED STUFF ========
rng ('default');
tt      = (0:tmax-1)' * timestep;
Xng     = rand ([nmax gmax]);     % initial conditions   % [.4    ; .6];
clc; close all;

%======== COMPUTE TRAJECTORIES ========
xntg = computeTrajsSimple (grn, Xng, tt);

%======== SET OPTIONS AND RUN infer ========
% NOTE: slopethresh, exprthresh, splinesmoothing are in pvxOpts_ngo
% NOTE: Rld_tsafety also should be removed.
pvxOpts_ngo = NaN (nmax, gmax, 3);
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
[grnFIGR, diagnostics] = infer (opts, xntg, tt, gmax);

%======== EXTRACT GENE STATES yk AND VELOCITIES vk ========
gTarget = 1;
kmax = nmax * tmax;
xkg = reshape (xntg, [kmax gmax]);
%ykg = reshape (diagnostics.yntg, [kmax gmax]);
%vkg = reshape (diagnostics.vntg, [kmax gmax]);
yk = diagnostics.yntg(:,:,gTarget); yk = yk (:);
vk = diagnostics.vntg(:,:,gTarget); vk = vk (:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======== VISUALIZE GENE STATES ========
close all;
set (groot, 'DefaultAxesFontSize', 16);

figure ('Position', [0 0 600 600], 'MenuBar', 'none'); hold on;
set (gcf, 'name', 'Data');
title ('Gene velocity v_A (green=positive, red=negative)'); %, 'interpreter','latex');
xlabel ('x_A');
ylabel ('x_B');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
%plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.4 .4 .4]); %gray
%======== Datapoints (filled circles) ========
gTarget = 1;
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




%======== COMPUTE MODEL PREDICTIONS FOR VELOCITIES (fitted v's) ========
% Use the MODEL (grnFIGR) to PREDICT the values of v at the ORIGINAL xkg points.
% 
% IDEALLY, WE WOULD JUST CALL velocityFunc
% BUT IT IS TOO DEEPLY NESTED
%
Tgg = grnFIGR.Tgg;
hg = grnFIGR.hg;
Rg = grnFIGR.Rg;
lambdag = grnFIGR.lambdag;
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
figure ('Position', [600 0 600 600], 'MenuBar', 'none'); hold on;
set (gcf, 'name', 'Fitted Model Predictions');
title ('Gene velocity v_A (green=positive, red=negative)'); %, 'interpreter','latex');
xlabel ('x_A');
ylabel ('x_B');
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




%==================================================================
%  computeTrajsSimple
%  Given initial conditions x_{ng}(t=0) = Xng,
%  compute trajectories xntg by solving diffusionless Glass equations
%
%    dx_{ng}/dt = v_{ng}  ({  x_{ng} }).
%==================================================================
function xntg = computeTrajsSimple (grn,Xng,tt)
nmax = size (Xng,1);
gmax = size (Xng,2);
tmax = numel (tt);
Xj = Xng(:);                             % pack initial conditions
[~,xtj] = ode45 (@(t,xj) velocityFunc(t,xj,grn), tt, Xj);  % ode45
xtng = reshape (xtj, [tmax nmax gmax]);  % unpack trajectories
xntg = permute (xtng, [2 1 3]);          % into correct format

%==================================================================
%  velocityFunc
%  Computes the Glass model gene velocity function, i.e.,
%  the RHS of the Glass ODE above, which is
%
%    v_{ng} = R_g S( sum_f T_gf x_nf + h_g ) - lambda_g x_ng.
%==================================================================
    function ddt_xj = velocityFunc (t,xj,grn)
        Rg = grn.Rg;    % Unpack grn structure here.
        Tgg = grn.Tgg;  % This makes the code easier to read,
        hg = grn.hg;    % but probably results in needless copying.
        lambdag = grn.lambdag;  % Note that nmax and gmax are inherited!
        
        xng = reshape (xj, [nmax gmax]);
        ddt_xng = (Rg .* heaviside (Tgg * xng' + hg))' - lambdag' .* xng;
        ddt_xj = reshape (ddt_xng, [nmax*gmax 1]);
    end
end
