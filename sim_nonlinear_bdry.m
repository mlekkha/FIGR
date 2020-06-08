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


%======== Define global structs for options and ODE options 
global opts;
global ODEopts;
global optimopts;
                        
%======== SET OPTIONS (see README.md for description) ========
opts = struct(  'debug', 0, ...                 
                'slopethresh', 0.01, ...         
                'exprthresh', 0.2, ...
                'splinesmoothing', 1.0, ...
                'spatialsmoothing', 0.5, ...
                'minborder_expr_ratio', 0.01, ...
                'Rld_method', 'slope_nodiff', ...
                'Rld_tsafety', 3, ...
                'synthesisfunction', 'synthesis_heaviside', ...
                'ODEAbsTol', 1e-4, ...
                'ODEsolver', 'ode45');

%======== set up ODE options 
ODEopts = odeset('AbsTol', opts.ODEAbsTol);

%======== SETUP GENE REGULATORY NETWORK PARAMETERS IN grn STRUCT ========
grn.Tgg     = [+0.0  0  +1 ; +0.5  +0.0  -1; 0  -1  +0.0];     % Regulatory coefficients
grn.hg      = [-0.23       0.2      0.55]';        % Thresholds/biases
grn.Rg      = [1           1        1]';           % Maximum synthesis rates
grn.lambdag = [4           4        4]';           % Degradation rates
grn.Dg      = [0           0        0]';          % Diffusion coefficients
xmax        = grn.Rg./grn.lambdag;                % maximum possible
                                                  %expression

%======== SETUP TIMEPOINTS tt ========
tt = [0 : 0.05: 3];  % Time points 

%======== SETUP INITIAL CONDITIONS xntg(:,1,:) ========
numNuclei     = 1;      % number of nuclei and/or experimental conditions
numTimepoints = 61;     % number of timepoints (= numel (tt))
numGenes      = 3;      % number of genes
xntg          = NaN (numNuclei,numTimepoints,numGenes); % allocate array
xntg(:,1,:)   = [.00  .00   0.05];                             % initial conditions

%======== COMPUTE TRAJECTORIES xntg ========
[xntg] = computeTrajs (opts, grn, xntg, tt);

%======== INFER GRN PARAMETERS grnFIGR to get yntg
[grnFIGR, diagnostics] = infer (opts, xntg, tt, numGenes);
yntg = diagnostics.yntg;

xntg_g2_on = xntg(:, yntg(:,:,2) > 0 , :);
xntg_g2_off = xntg(:, yntg(:,:,2) < 0 , :);

% compute the surface according to equations
x1 = [0: 0.01*xmax(1): xmax(1)]';
x2 = [0: 0.01*xmax(2): xmax(2)]';
x3 = -1.0*([x1 x2]*grn.Tgg(2,1:2)' + ones(size(x1,1),1)*grn.hg(2))/grn.Tgg(2,3);
Z = repmat(x3', size(x1,1) ,1);


%======== VISUALIZE RESULTS
close all;
set (gcf, 'Position', [0 0 800 400]);


subplot (2, 2, 1); 
%-------- PLOT TRAJECTORIES x1(t) AND x2(t) --------
plot (tt, squeeze(xntg(1,:,:))); xlabel ('t'); legend ('x_1', 'x_2', 'x_3');
title('Synthetic data');

%-------- PLOT TRAJECTORIES IN (x1,x2) PLANE --------
subplot (2, 2, 2); hold on;
scatter (xntg_g2_on(1,:,1), xntg_g2_on(1,:,2), 'g'); xlabel ('x_1'); ylabel ('x_2');
scatter (xntg_g2_off(1,:,1), xntg_g2_off(1,:,2), 'r'); xlabel ('x_1'); ylabel ('x_2');
plot(x1, x2, 'b');
%-------- PLOT TRAJECTORIES IN (x1,x3) PLANE --------
subplot (2, 2, 3); hold on;
scatter (xntg_g2_on(1,:,1), xntg_g2_on(1,:,3), 'g'); xlabel ('x_1'); ylabel ('x_3');
scatter (xntg_g2_off(1,:,1), xntg_g2_off(1,:,3), 'r'); xlabel ('x_1'); ylabel ('x_3');
plot(x1, x3, 'b');
%-------- PLOT TRAJECTORIES IN (x2,x3) PLANE --------
subplot (2, 2, 4); hold on;
scatter (xntg_g2_on(1,:,2), xntg_g2_on(1,:,3), 'g'); xlabel ('x_2'); ylabel ('x_3');
scatter (xntg_g2_off(1,:,2), xntg_g2_off(1,:,3), 'r'); xlabel ('x_2'); ylabel ('x_3');
plot(x2, x3, 'b');


%-------- PLOT THEORETICAL SWITCHING BOUNDARIES IN 3D --------
% for gene 2
figure;
hold on;
scatter3(xntg_g2_on(:,:,1), xntg_g2_on(:,:,2), xntg_g2_on(:,:,3), 'g');
scatter3(xntg_g2_off(:,:,1), xntg_g2_off(:,:,2), xntg_g2_off(:,:,3), 'r');
surf(x1, x2, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
view(0, 0.5729);

