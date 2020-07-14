%======================================
% MATLAB SCRIPT
% - Make "velocity plot"
%========================================

%======== USER-EDITABLE PARAMETERS ========
gmax     = 2;      % 2-gene toy model
fmax     = 2;      % EXTERNAL REGULATORS NOT SUPPORTED!
nmax     = 20;     % number of nuclei (initial coords)
tmax     = 11;     % number of timepoints
timestep = 0.1;

gTarget = 2;         % gene whose VELOCITY we want to visualize
colorON = [0 .5 0];  % green indicates v>0
colorOFF = [.7 0 0]; % red indicates v>0

rng ('shuffle');     % randomize.     OR  rng ('default');

tt      = (0:tmax-1)' * timestep;


%======== VERSION 1: PREDETERMINED GRN
% grn.Tgg     = [1 -.5 ; -.5 1];   % autoactivating and mutually repressing
% grn.hg      = [-.25  ; -.25];    % leads to tetrastable behavior
% grn.Rg      = [2     ; 1];
% grn.lambdag = [2     ; 1];
% 
% Xng     = rand ([nmax gmax]);     % initial conditions   % [.4    ; .6];

%======== VERSION 2: RANDOM GRN
%======== GENERATE RANDOM PARAMETERS T,h
%	Choose a point r0 within the unit hypercube using rand().
%	Choose a unit random vector T in N=(gmax+numExternals) dimensions,
%		by generating a vector of N normal deviates and normalizing it.
%	Consider the hyperplane T.(r - r0) = 0, that is, T.r + h = 0
%		where h = -T.r0.
grn.Tgg = NaN (gmax, fmax);
grn.hg = NaN (gmax, 1);
for g=1:gmax
    r0vec    = rand ([fmax 1]);
    Tvec     = normrnd (0, 1, [fmax 1]);
    Tvec     = Tvec ./ norm(Tvec);
    h        = -Tvec' * r0vec;
    grn.Tgg(g,:) = Tvec';
    grn.hg(g)    = h;
end   % The above could be modified so that T is auto-activating or auto-repressing
%======== GENERATE RANDOM PARAMETERS R,lambda
grn.Rg      = rand (gmax, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
grn.lambdag = rand (gmax, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
grn.Dg      = zeros (gmax, 1);
%======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
Xng         = rand (nmax, fmax) .* grn.Rg' ./ grn.lambdag';


%======== DERIVED STUFF ========
clc; close all;

%======== COMPUTE TRAJECTORIES ========
xntg = computeTrajsSimple (grn, Xng, tt);
% xntg        = NaN ([nmax,tmax,fmax]);
% xntg(:,1,:) = Xng;
% xntg(:,:,gmax+1:end) = 0.;
% xntg = computeTrajs (opts, grn, xntg, tt); 
%  %%% original was somethnig like this

%======== SET OPTIONS AND RUN infer ========
% NOTE: Rld_tsafety also should be removed.
pvxOpts_ngo = NaN (nmax, gmax, 3);
pvxOpts_ngo(:,:,1) = 1.00;  % p (spline unsmoothing parameters)
pvxOpts_ngo(:,:,2) = 0.01;  % v (velocity thresholds)
pvxOpts_ngo(:,:,3) = 0.20;  % x (expression thresholds)
opts = struct(  'debug', 0, ...
    'Rld_tsafety', 0, ...       % USUALLY 3
    'spatialsmoothing', 0.5, ...
    'minborder_expr_ratio', 0.01, ...
    'Rld_method', 'slope', ...
    'synthesisfunction', 'synthesis_sigmoid_sqrt', ...
    'ODEAbsTol', 1e-5, ...      % originally 1e-3
    'ODEsolver', 'ode45', ...
    'pvxOpts_ngo', pvxOpts_ngo, ...
    'lm', 'FIGRlogReg', 'lambda', 0.001);  % Joanna used lambda = 0.5,

%'lm', 'glmfit'); % glmfit|FIGRlogReg|lassoglm

% 
% FIGRlogReg is broken probably becasue the derivatives are wrong

[grnFIGR, diagnostics] = infer (opts, xntg, tt, gmax);


%======== EXTRACT GENE STATES yk AND VELOCITIES vk ========
kmax = nmax * tmax;
xkg = reshape (xntg, [kmax gmax]);
ykg = reshape (diagnostics.yntg, [kmax gmax]);
vkg = reshape (diagnostics.vntg, [kmax gmax]);
yk = ykg(:,gTarget);
vk = vkg(:,gTarget);
vkDATA = vk; % save for later use

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======== VISUALIZE GENE STATES ========
close all;
set (groot, 'DefaultAxesFontSize', 16);

figure ('Position', [0 0 400 400], 'MenuBar', 'none'); hold on;
set (gcf, 'name', 'Data');
title ('Gene velocity v_B^{data}'); %, 'interpreter','latex');
xlabel ('x_A');
ylabel ('x_B');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.7 .7 .7]); %gray
%======== Datapoints (filled circles) ========
markerSizes  = 100.0 * abs (vk) .^ 0.5;  % bigger markers mean stronger vels
markerColors = NaN (kmax, 3);
markerColors = zeros (kmax, 3);
kSel = find(yk>0);
markerColors(kSel,:) = markerColors(kSel,:)*0 + colorON;
scatter (xkg(kSel,1), xkg(kSel,2), markerSizes(kSel), markerColors(kSel,:), '^','filled');
kSel = find(yk<0);
markerColors(kSel,:) = markerColors(kSel,:)*0 + colorOFF;
scatter (xkg(kSel,1), xkg(kSel,2), markerSizes(kSel), markerColors(kSel,:), 'v');
%======== Switching boundary ========
for g=gTarget:gTarget
    f = @(x1,x2)   grn.Tgg(g,1)*x1 + grn.Tgg(g,2)*x2 + grn.hg(g);
    fcontour (f, [0 1 0 1], 'LevelList', [0], 'LineStyle', '--', 'LineColor', 'b', 'LineWidth', 1);
end
%======== Legend ========
handles(1) = scatter(NaN,NaN,9,colorON,'^','filled');
handles(2) = scatter(NaN,NaN,9,colorOFF,'v');
legend (handles, 'v_B > 0', 'v_B < 0'); legend 



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
vkg = (Rg .* heaviside (Tgg * xkg' + hg))' - lambdag' .* xkg;
ykg = sign (vkg);   
yk = ykg(:,gTarget);
vk = vkg(:,gTarget);
vkMODEL = vk; % save for later use



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure ('Position', [400 0 400 400], 'MenuBar', 'none'); hold on;
set (gcf, 'name', 'Fitted Model Predictions');
title ('Gene velocity v_B^{model}'); %, 'interpreter','latex');
xlabel ('x_A');
ylabel ('x_B');
axis ([0 1 0 1]);
%======== Trajectories (curves) ========
plot (squeeze (xntg(:,:,1))', squeeze (xntg(:,:,2))', 'color', [.7 .7 .7]); %gray
%======== Datapoints (filled circles) ========
markerSizes  = 100.0 * abs (vk) .^ 0.5;  % bigger markers mean stronger vels
markerColors = zeros (kmax, 3);
kSel = find(yk>0);
markerColors(kSel,:) = markerColors(kSel,:)*0 + colorON;
scatter (xkg(kSel,1), xkg(kSel,2), markerSizes(kSel), markerColors(kSel,:), '^','filled');
kSel = find(yk<0);
markerColors(kSel,:) = markerColors(kSel,:)*0 + colorOFF;
scatter (xkg(kSel,1), xkg(kSel,2), markerSizes(kSel), markerColors(kSel,:), 'v');
%======== Switching boundary (inferred) ========
for g=gTarget:gTarget
    f = @(x1,x2)   grnFIGR.Tgg(g,1)*x1 + grnFIGR.Tgg(g,2)*x2+ grnFIGR.hg(g);
    fcontour (f, [0 1 0 1], 'LevelList', [0], 'LineStyle', '--', 'LineColor', 'b', 'LineWidth', 1);
end
%======== Legend ========
handles(1) = scatter(NaN,NaN,9,colorON,'^','filled');
handles(2) = scatter(NaN,NaN,9,colorOFF,'v');
legend (handles, 'v_B > 0', 'v_B < 0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:gmax
    mynorm          = norm( grn.Tgg(g,:) );   % normalize T_original
    grn.Tgg(g,:)    = grn.Tgg(g,:) / mynorm;
    grn.hg(g)       = grn.hg(g)    / mynorm;
    mynorm          = norm( grnFIGR.Tgg(g,:) );% normalize T_inferred
    grnFIGR.Tgg(g,:)    = grnFIGR.Tgg(g,:) / mynorm;
    grnFIGR.hg(g)       = grnFIGR.hg(g)    / mynorm;
end

figure ('Position', [0 400 400 400], 'MenuBar', 'none'); hold on;
axis (10*[-1 1 -1 1]); grid ON;
title ('T-vector misfit plot');
xlabel ('T_{B<-A,B}^{original}');
ylabel ('T_{B<-A,B}^{inferred}');
scatter (grn.Tgg(gTarget,:), grnFIGR.Tgg(gTarget,:), 'filled');
plot (999*[-1 1], 999*[-1 1], 'color', [.4 .4 .4]);  % "y=x" line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure ('Position', [400 400 400 400], 'MenuBar', 'none'); hold on;
axis ([-1 1 -1 1]); grid ON;
title ('Velocity misfit plot');
xlabel ('v_B^{data}');
ylabel ('v_B^{model}');
scatter (vkDATA, vkMODEL);
plot ([-1 1], [-1 1], 'color', [.4 .4 .4]);  % "y=x" line
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
