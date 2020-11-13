function fh_roc = neyda(opts, grn, ExtInpInterp)

% Returns handle to nested function that computes the rate of change
% according to the gene circuit equations
%
%    d/dt x_ntg = R_g g( sum_f T_gf x_ntf + h_g) - lambda_g x_ntg
%                       + D_g(x_(n+1)tg + x_(n-1)tg - 2 x_ntg)
%

numGenes = numel(grn.Rg);                  % G
numRegulators = size(grn.Tgg, 2);          % F
numExternals = numRegulators - numGenes;   % F-G

% local variables for the nested function.
% MODIFIED BY YLL 2018-12-11

Tgf = grn.Tgg;   % INCLUDE ALL GENES
hg = grn.hg;
Rg = grn.Rg;
lambdag = grn.lambdag;
Dg = grn.Dg;

if (numExternals > 0)
    tt         = ExtInpInterp.tt;
    xk_ext_dat = ExtInpInterp.xk_ext_dat;
else
    tt = 0;
    xk_ext_dat = 0;
end    


fh_roc = @getRateofChange;




%======================================================================
%             function vk = getRateofChange(t, xk)
%
% Nested function that computes gene expression "velocities" vk 
% given "positions" xk (gene product concentrations) and current time t.
% xk is a column vector, flattened from xgn (a GxN matrix), 
% as required by MATLAB's ODE solver functions.
%
% The diffsion flux of protein g from nucleus n TO nucleus n+1 is
%
%   q_gn = D_g * (x_gn - x_{g,n+1}) .
%
% The rate of change of concentration of protein g in nucleus n is
%
%   d/dt (x_ng)
%   = v_ng
%   =   R_g fsigmoid (sum_f T_gf x_ntf + h_g)       (synthesis)
%     + q_{g,n-1} - q_gn                            (diffusion)
%     - lambda_g x_ntg                              (degradation)
%
% where 
%   g=1,...,G      G = numGenes      = 4 (HKGN)
%   f=1,...,F      F = numRegulators = 7 (HKGNBCT)
%   n=1,...,N      N = numNuclei     = 58
%
% INPUTS:
%   t,xk
% VARIABLES INHERITED FROM ENCLOSING SCOPE:
%   numGenes
%   Tgg,hg,Rg,lambdag,Dg
%   tt,xk_ext_dat
% OUTPUTS: 
%   VK
%
% YLL 2018-12-11: Reshaping to 2D matrix for readability.

    function vk = getRateofChange(t, xk)   % THE LHS IS "dx/dt"
        
        numEqns = numel(xk);
        numNuclei = numEqns/numGenes;

        % initialize arrays

        %---- "INTERNAL" GENE CONCENTRATIONS xgn
        xgn = reshape (xk, [numGenes numNuclei]);

        % concentrations for interior product with Tgf
        xfn = nan (numGenes+numExternals, numNuclei);

        % accumulates rate of change
        vgn = nan (numGenes, numNuclei);

        % flattened rate of change to return
        vk = nan (numEqns, 1);
 
        % flux from diffusion
        qgn = nan (numGenes, numNuclei);

        % external input concentrations
        wgn = nan (numExternals, numNuclei);

        % flattened external input concentrations
        xk_ext = nan (numExternals*numNuclei, 1);

        %---- "EXTERNAL" GENE CONCENTRATIONS wgn
        if (numExternals > 0)
            xk_ext = interp1(tt, xk_ext_dat', t)';
            wgn = reshape (xk_ext, [numExternals numNuclei]);

            % include both internal and external concentrations for Tgf*xf
            xfn = [xgn; wgn];      
        else
            wgn = zeros(0, numNuclei);

            % only internal concentrations for Tgf*xf
            xfn = xgn;      
        end

        % Flux Dg*(x_(n) - x_(n+1))
        % NOT using implicit expansion since MATLAB Coder does not support
        % it
        qgn = repmat(Dg,1,numNuclei) .* (xgn - circshift (xgn, numNuclei-1, 2));  % rotate LEFT in "n"

        qgn(:,numNuclei) = zeros(numGenes, 1); % NO FLUX FROM LAST NUCLEUS TO FIRST ONE
        
        % Synthesis
        % Since arrayfun not supported by MATLAB Coder,
        % evaluate synthesis in each nucleus separately
        for j=1:numNuclei

            for k=1:numGenes

                % which synthesis function we should use
                % can't use variable function handles since MATLAB code 
                % generation does support it
                if strcmp(opts.synthesisfunction, 'synthesis_sigmoid_sqrt')
                    vgn(k,j) = Rg(k) * synthesis_sigmoid_sqrt(Tgf(k,:) * xfn(:,j) + hg(k));
                elseif strcmp(opts.synthesisfunction, 'synthesis_heaviside')
                    vgn(k,j) = Rg(k) * synthesis_heaviside(Tgf(k,:) * xfn(:,j) + hg(k));
                else
                    error('Unknown synthesis function %s',...
                                opts.synthesisfunction);
                end
   
            end
            
        end


        % Diffusion
        % second term is Dg*(x_(n-1) - x_(n)), third is Dg*(x_(n) - x_(n+1))
        vgn = vgn + circshift (qgn, 1, 2) - qgn;            

        % Degradation
        % NOT using implicit expansion since MATLAB Coder does not support
        % it
        vgn = vgn - repmat(lambdag,1,numNuclei) .* xgn;             

        % flatten to return
        vk = reshape (vgn, numEqns, 1);
        
    end % from getRateofChange()

end % from RateofChange()

function g = synthesis_sigmoid_sqrt(u)
g = 0.5*(u/sqrt(1.0+u*u) + 1.0);
end % from synthesis_sigmoid_sqrt()

function g = synthesis_heaviside(u)
if u >= 0.0 % Heaviside function,
    g = 1.0; %differs from inbuilt MATLAB one in that g(u) = 1
else         % instead of 0.5 if u = 0.
    g = 0.0;
end
end % from synthesis_heaviside()
