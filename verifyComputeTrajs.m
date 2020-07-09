%======================================
% MATLAB SCRIPT
% - verifies that computeTrajs works
% - provides computeTrajsSimple (neat, but not as general)
%========================================

%========================================
% NOTES ON INDEX PACKING
%
% For the purposes of FIGR inference via logistic regression,
% it is convenient to define the DATAPOINT index k=[n,t]
% as a composite of the NUCLEUS index n and the TIME index t.
%
% For the purposes of computeTrajs,
% we define the DYNAMICAL VARIABLE index j=[n,g]
% as a composite of the NUCLEUS index n and the GENE index g,
% because the "dynamical variable" in ode45 must be a column vector.
%
% kmax   = nmax * tmax;
% jmax   = nmax * gmax;
%========================================

%======== USER-EDITABLE PARAMETERS ========
gmax     = 2;      % 2-gene toy model
fmax     = 2;      % EXTERNAL REGULATORS NOT SUPPORTED!
nmax     = 20;     % number of nuclei (initial coords)
tmax     = 11;     % number of timepoints
timestep = 0.1;
grn.Tgg     = [1 -.5 ; -.5 1];   % autoactivating and mutually repressing
grn.hg      = [-.25  ; -.25];    % leads to tetrastable behavior
grn.Rg      = [2     ; 1];
grn.lambdag = [2     ; 1];

%======== DERIVED STUFF ========
rng ('default');
tt      = (0:tmax-1)' * timestep;
Xng     = rand ([nmax gmax]);     % initial conditions   % [.4    ; .6];
clc; close all;

%======== COMPUTE TRAJECTORIES ========
xntg = computeTrajsSimple (grn, Xng, tt);      % CALL LOCAL computeTrajsSimple

% grn.Dg = [0 ; 0];
% opts.debug = 0;
% opts.ODEsolver = 'ode45';
% opts.synthesisfunction = 'synthesis_heaviside';
% %opts.synthesisfunction = 'synthesis_sigmoid_sqrt';
% Xntg    = NaN (nmax,tmax,gmax); Xntg(:,1,:) = Xng;
% xntg = computeTrajs (opts, grn, Xntg, tt);  % CALL ORIGINAL computeTrajs


%======== FIGURE 1: TRAJECTORIES IN (x_1, x_2) PLANE ========
figure ('Position', [0 0 600 600]); hold on;
title ('TRAJECTORIES IN  (x_1, x_2) PLANE');
xlabel ('x_A'); ylabel ('y_A'); axis([0 1 0 1]);
scatter (Xng(:,1), Xng(:,2), 'ro');         % plot initial conds
for n=1:nmax
    plot (xntg(n,:,1), xntg(n,:,2), 'bx-'); % plot trajs
end

%======== FIGURE 2: TRAJECTORIES x_2(t) ========
figure ('Position', [600 0 600 600]); hold on;
title ('FIGURE 2: TRAJECTORIES x_2(t)');
scatter (0*Xng(:,2), Xng(:,2), 'ro');         % plot initial conds
for n=1:nmax
    plot (tt, xntg(n,:,2), 'bx-');
end

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
