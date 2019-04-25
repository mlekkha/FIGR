function [xntg] = computeTrajs (opts,grn,Xntg,tt)

% This function solves the differential equations
%
%   d/dt x_ntg = R_g g( sum_f T_gf x_ntf + h_g) - lambda_g x_ntg
%                       + D_g(x_(n+1)tg + x_(n-1)tg - 2 x_ntg)
%
% where x_ntg for n=1:numGenes                        are provided only at t=0
% and   x_ntg for n=numGenes+1:numGenes+numExternals  are provided at all t values.
%
% ARGUMENTS:
%   grn.Tgg            genetic interconnect matrix
%   grn.hg             thresholds
%   grn.Rg             maximum synthesis rates
%   grn.lambdag        degradation rates
%   grn.Dg             diffusion constants
%   Xntg(:,1,1:F)      initial conditions for all regulators (int & ext)
%   Xntg(:,:,G+1:F)    trajectories for all external regulators
%   tt                 timepoints at which to store the trajectories
% RETURNS:
%   xntg(:,:,:)        trajectories for all genes
%
% NOTES:
%   See Example Program 2 for an illustration of Xntg and xntg 
%    when upstream regulators are present.
%

if opts.debug > 0
    fprintf ('======== computeTrajs() ========\n');
    fprintf ('======== VALIDATING INPUTS ========\n');
end

numNuclei     = size (Xntg, 1);
numTimepoints = size (Xntg, 2);
numGenes      = numel (grn.Rg);
numExternals  = size (Xntg,3) - numGenes;

if opts.debug > 0
    fprintf ('tt = \n %s\n'    , mat2str (tt') );
    fprintf ('numNuclei          = %d \n',numNuclei);
    fprintf ('numTimepoints      = %d \n',numTimepoints);
    fprintf ('numGenes           = %d \n',numGenes);
    fprintf ('numExternals       = %d \n',numExternals);
    fprintf ('size(Xntg)         = %d %d %d \n', size(Xntg));
end

assert ( isequal ( size(grn.Tgg    ), [numGenes   numGenes+numExternals] ) );
assert ( isequal ( size(grn.hg     ), [numGenes, 1] ) );
assert ( isequal ( size(grn.Rg     ), [numGenes, 1] ) );
assert ( isequal ( size(grn.lambdag), [numGenes, 1] ) );
assert ( isequal ( size(grn.Dg     ), [numGenes, 1] ) );

% if we have external inputs, then flatten Xntg so that it may be
% passed to RateOfChange().
%
% YLL: EVENTUALLY I MAY WISH TO MAKE xk_ext_dat JUST wgn_dat OR SOMETHING
% LIKE THAT
if (numExternals > 0)
    
    ExtInpInterp.tt = tt;
    
    % Shift dimensions so that genes change with rows, nuclei with
    % columns, and time is the third dimension
    Wgnt = shiftdim(Xntg(:,:,numGenes+1:numGenes+numExternals), 2);
    
    % Flatten
    ExtInpInterp.xk_ext_dat = reshape(Wgnt, ...
        [numExternals*numNuclei, numTimepoints]);
    
    % Plot patterns to ensure that the flattening hasn't clobbered Xntg
    if (opts.debug > 2)
        figure;
        for j=1:numTimepoints
            
            hold on;
            for k=1:numExternals
                plot(1:numNuclei, ...
                    ExtInpInterp.xk_ext_dat(k:numExternals:end,j));
            end
            pause;
            hold off; clf;            
        end
    end
    
else
    
    ExtInpInterp.tt = 0;
    ExtInpInterp.xk_ext_dat = 0;
    
end % from if (numExternals > 0)

% Flatten Xntg so we can provide initial conditions
Xgnt = shiftdim(Xntg(:,:,1:numGenes), 2);
Xkt = reshape(Xgnt, numGenes*numNuclei, numTimepoints);

% Initialize the RateOfChange function and save the "odefun" style
% function that is returned.
roc = RateofChange(opts, grn, ExtInpInterp);

% set ODE opts
odeOpts = odeset('AbsTol', opts.ODEAbsTol);

% solverfunc
solver = str2func(opts.ODEsolver);

% solve; note that in the returned solution, time varies with row not
% with column as in Xkt
[tts, xtk] = solver(roc, tt, Xkt(:,1), odeOpts);

% reshape solution to return xntg
xtgn = reshape(xtk, [numel(tts), numGenes, numNuclei]);

% YLL 2018-12-11
% PROBLEM: 
% According to the MATLAB documentation for reshape(),
% Beyond the second dimension, the output, B, does not reflect trailing dimensions with a size of 1. 
% For example, reshape(A,3,2,1,1) produces a 3-by-2 matrix.
% WORKAROUND: INSTEAD OF WRITING
%      xntg = shiftdim(xtgn, 2);
%   CODE THIS BY HAND AS FOLLOWS (HORRIBLY INEFFICIENT):
for n=1:numNuclei
    for t=1:numTimepoints
        for g=1:numGenes
            xntg(n,t,g) = xtgn(t,g,n);
        end
    end
end
%disp (size(xntg)); 

% YLL 2018-11-19
%======= xntg contains trajectories for HKGN.  Now tack on the BCT trajs.
% In general we may need to interpolate the BCT trajs to the desired timepoints...
xntg = cat (3, xntg, Xntg(:,:,numGenes+1:end));

end % from computeTrajs()
