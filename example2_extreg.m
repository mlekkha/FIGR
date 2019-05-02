clc;
fprintf (['\n'...
    '============================== Example 2 ==============================\n'...
    'This example studies a toy model with 2 genes plus 1 upstream regulator.\n'...
    'The expression of the upstream regulator increases with time from 0.2 to 1.7.\n'...
    'Genes 1 and 2 have almost no influence on each other.\n'...
    'Gene 3 activates gene 1 and represses gene 2.\n'...
    'Thresholds are such that gene 1 switches ON when x_3 > 0.6 \n'...
    ' and gene 2 switches OFF when x_3 > 1.2. \n'...
    'In this case, there is only one trajectory, and only two switching points, \n'...
    '  so FIGR is unable to infer correct values for regulatory parameters (T,h); \n'...
    '  moreover, FIGR can only infer T and h up to some multiplicative constant. \n'...
    'Nevertheless, kinetic parameters (R,lambda) are inferred well, \n'...
    '  and in this case the inferred GRN reproduces the original trajectories fairly well. \n'...
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

%======== Define global struct for ODE options and set it up
global ODEopts;
ODEopts = odeset('AbsTol', opts.ODEAbsTol);

%======== SETUP GENE REGULATORY NETWORK PARAMETERS IN grn STRUCT ========
grn.Tgg     = ...
    [0      +0.01   -1  ;...
    -0.01    0      +1  ];        % Regulatory coefficients (2x3 matrix)

grn.hg      = [+0.6      -1.2]';        % Thresholds/biases
grn.Rg      = [1           1]';           % Maximum synthesis rates
grn.lambdag = [1           1]';           % Degradation rates
grn.Dg      = [0           0 ]';          % Diffusion coefficients

%======== SETUP TIMEPOINTS tt ========
tt = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1. 1.1 1.2 1.3 1.4 1.5]';  % Time points

%======== SETUP INITIAL CONDITIONS xntg(:, 1, 1:2)
%======== AND EXTERNAL REGULATORS' TRAJECTORIES xntg(:, :, 3)
numNuclei     = 1;           % number of nuclei and/or experimental conditions
numTimepoints = numel(tt);   % number of timepoints
numGenes      = 2;                     % number of "internal" genes in network
numExternals  = 1;                     % number of upstream regulators
numRegulators = numGenes+numExternals; % total number of genes
xntg = NaN(numNuclei, numTimepoints, numRegulators);
xntg(1,:,:) = [...
   .96  .03 .2 ; ...  % initial conditions (first timepoint)
    0   0   .3 ; ...
    0   0   .4 ; ...
    0   0   .5 ; ...  % need to provide time series data for
    0   0   .6 ; ...  % upstream regulator
    0   0   .7 ; ...
    0   0   .8 ; ...
    0   0   .9 ; ...  % the "0" elements in this matrix
    0   0   1 ; ...   % are not needed by computeTrajs()
    0   0   1.1 ; ...
    0   0   1.2 ; ...
    0   0   1.3 ; ...
    0   0   1.4 ; ...
    0   0   1.5 ; ...
    0   0   1.6 ; ...
    0   0   1.7];

%======== COMPUTE TRAJECTORIES xntg ========
% This overwrites the "0" placeholders in xntg.
[xntg] = computeTrajs (opts, grn, xntg, tt);

%======== INFER GRN PARAMETERS grnFIGR USING FIGR ALGORITHM
[grnFIGR, diagnostics] = infer (opts, xntg, tt, numGenes);

%======== RESIMULATE TRAJECTORIES xntgFIGR BASED ON INFERRED GRN PARS
[xntgFIGR] = computeTrajs (opts, grnFIGR, xntg, tt);


%======== DISPLAY grnSA AND grnFIGR FOR COMPARISON
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

%-------- PLOT ORIGINAL AND RESIMULATED TRAJECTORIES IN (x1,x2) PLANE --------
subplot (1, 2, 1); hold on;
plot (xntg    (1,:,1), xntg(1,:,2), 'ro-');
plot (xntgFIGR(1,:,1), xntgFIGR(1,:,2), 'ro--');
xlabel ('x_1'); ylabel ('x_2'); axis ([0 1 0 1]);

%-------- PLOT ORIGINAL AND RESIMULATED TRAJECTORIES x1(t), x2(t), x3(t) --------
subplot (1, 2, 2); hold on;
colorOrder = get (gca,'colororder');
for g=1:3
    plot (tt, squeeze(xntg    (1,:,g)), 'x-',  'Color', colorOrder(g,:));
end
for g=1:3
    plot (tt, squeeze(xntgFIGR(1,:,g)), 'x--', 'Color', colorOrder(g,:));
end
xlabel ('t'); legend ('x_1', 'x_2', 'x_3');
