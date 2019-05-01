clc;
fprintf (['\n'...
'============= Example 1 (Figure 1 of Fehr et al.  manuscript) =============\n'...
'In this example, a Gene Regulatory Network (GRN) is defined by parameters T,h,R,lambda,D. \n'...
'Starting from one initial configuration, the gene expression trajectory x_g(t) is calculated. \n'...
'From this trajectory, FIGR is used to infer the parameters of the GRN. \n'...
'The output graphs show the trajectory (red circles), \n'...
'  actual switching boundaries (solid blue lines), \n'...
'  and inferred switching boundaries (dashed blue lines). \n'...
'This indicates that FIGR is able to infer parameters reasonably well even from \n'...
'  the small amount of data here. \n'...
'=======================================================================\n'...
'\n']);

%======== SET OPTIONS FOR COMPUTETRAJS ========
opts.synthesisfunction = 'synthesis_heaviside'; % simulate Glass model in Eq. (3)
opts.ODEAbsTol = 1e-4;                          % ODE solver tolerance
opts.ODEsolver = 'ode45';                       % 4th order Runge-Kutta
opts.debug = 0;                                 % verbosity level (0-3)
%======== SET OPTIONS FOR FIGR ========
opts.slopethresh = .01;            % velocity threshold v^c in Eq. (10)
opts.exprthresh = .2;              % expression threshold x^c in Eq. (10)
opts.splinesmoothing = 1.00;       % 1=none, 0=extreme smoothing
opts.Rld_method = 'slope_nodiff';  % inference method for kinetic parameters 
opts.Rld_tsafety = 3;              % points to ignore near switching events


%======== SETUP GENE REGULATORY NETWORK PARAMETERS IN grn STRUCT ========
grn.Tgg     = [-.1   +1 ;  -1.    0];     % Regulatory coefficients
grn.hg      = [-0.45       0.35]';        % Thresholds/biases
grn.Rg      = [1           2]';           % Maximum synthesis rates
grn.lambdag = [1           2]';           % Degradation rates
grn.Dg      = [0           0 ]';          % Diffusion coefficients

%======== SETUP TIMEPOINTS tt ========
tt = [0 : 0.05: 3];  % Time points 

%======== SETUP INITIAL CONDITIONS xntg(:,1,:) ========
numNuclei     = 1;      % number of nuclei and/or experimental conditions
numTimepoints = 61;     % number of timepoints (= numel (tt))
numGenes      = 2;      % number of genes
xntg          = NaN (numNuclei,numTimepoints,numGenes); % allocate array
xntg(:,1,:)   = [.05  .00];                             % initial conditions

%======== COMPUTE TRAJECTORIES xntg ========
[xntg] = computeTrajs (opts, grn, xntg, tt);

%======== INFER GRN PARAMETERS grnFIGR USING FIGR ALGORITHM
[grnFIGR, diagnostics] = infer (opts, xntg, tt, numGenes);
%yntg = diagnostics.yntg;


%======== COMPUTE TRAJECTORIES xntgFIGR from inferred  GRN PARAMETERS 
[xntgFIGR] = computeTrajs (opts, grnFIGR, xntg, tt);

%======== DISPLAY grn AND grnFIGR FOR COMPARISON
% Since Tgg and hg can only be determined to within a multiplicative
% factor, normalize so that the length of the Tg vector is 1.

% Normalize actual GRN grn
for g=1:numGenes
    fac = norm(grn.Tgg(g,:));
    grn.Tgg(g,:) = grn.Tgg(g,:) / fac;
    grn.hg(g,:)  = grn.hg(g,:) / fac;
end
disp ('GRN parameters of toy model: ');
disp (struct2table (grn));
disp ('');

% Normalize inferred GRN grnFIGR
for g=1:numGenes
    fac = norm(grnFIGR.Tgg(g,:));
    grnFIGR.Tgg(g,:) = grnFIGR.Tgg(g,:) / fac;
    grnFIGR.hg(g,:)  = grnFIGR.hg(g,:) / fac;
end
disp ('GRN parameters inferred using FIGR: ');
disp (struct2table (grnFIGR));


%======== VISUALIZE RESULTS
close all;
set (gcf, 'Position', [0 0 800 400]);


subplot (1, 2, 1); 
%-------- PLOT TRAJECTORIES x1(t) AND x2(t) --------
plot (tt, squeeze(xntg(1,:,:)), 'x' ); xlabel ('t'); legend ('x_1', 'x_2');
title('Synthetic data');

subplot (1, 2, 2); hold on;
%-------- PLOT TRAJECTORIES IN (x1,x2) PLANE --------
plot (xntg(1,:,1), xntg(1,:,2), 'r-'); xlabel ('x_1'); ylabel ('x_2');
%-------- PLOT THEORETICAL AND INFERRED SWITCHING BOUNDARIES --------
for g=1:2
    f = @(x,y)   grn.Tgg(g,1)*x + grn.Tgg(g,2)*y + grn.hg(g);
    fcontour (f, [0 1 0 1], 'LevelList', [0], 'LineColor', 'b', 'LineWidth', 1);
    f = @(x,y)   grnFIGR.Tgg(g,1)*x + grnFIGR.Tgg(g,2)*y + grnFIGR.hg(g);
    fcontour (f, [0 1 0 1], 'LevelList', [0], 'LineStyle', '--', 'LineColor', 'b', 'LineWidth', 1);
end
legend('Synthetic data trajectory', 'Actual switching hyperplane', 'Inferred switching hyperplane');
title('Inference');

