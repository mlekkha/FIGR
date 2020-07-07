function [grn, diagnostics] = infer (opts, xntg,tt,numGenes)
%
% This function implements the FIGR algorithm:
% 1. On/off states are guessed for each gene, in each nucleus, at each timepoint.
% 2. Regulatory parameters {T,h} are inferred using logistic regression.
% 3. Kinetic parameters {R,lambda,D} are inferred using linear regression
%      ('slope' method) or or methods as described in the paper.
%
%
% ARGUMENTS:
%   xntg                        gene expression trajectories x(n,t,g) as 3D array
%   tt                          timepoints as column vector
%   opts.slopethresh            velocity threshold v^c in Eq. (10)
%   opts.exprthresh             expression threshold x^c in Eq. (10)
%   opts.splinesmoothing        1=no smoothing, 0=extreme smoothing
%   opts.Rld_method             'slope', 'slope_nodiff', 'conc', or 'kink'
%   opts.Rld_tsafety            PARAMETER USED BY SLOPE METHOD
%   opts.minborder_expr_ratio   PARAMETER USED BY KINK METHOD
%   opts.spatialsmoothing       PARAMETER USED BY KINK METHOD
%   opts.geneNames              cell array (unimportant; mainly for plotting)
% RETURN VALUES:
%   grn.Tgg            genetic interconnect matrix
%   grn.hg             thresholds
%   grn.Rg             maximum synthesis rates
%   grn.lambdag        degradation rates
%   grn.Dg             diffusion constants
%   diagnostics.yntg   on/off states y(n,t,g)
%
% NOTES:
%   Indices generally run from n=1:N, t=1:T, g=1:G where
%   N=number of nuclei, T=number of timepoints, and G=number of genes.  
%   N, T, and G are determined from the dimensions of the "xntg" array.
%
%   For the Drosophila dataset, xntg contains data for 7 genes.
%   The infer routine is called with numGenes=4 to indicate that
%     we only wish to infer regulatory parameters acting on the first 4
%     genes; the last 3 genes are upstream regulators.

numNuclei     = size (xntg,1);
numTimepoints = size (xntg,2);
numRegulators = size (xntg,3);
numExternals = numRegulators - numGenes;
numDatapoints = numNuclei * numTimepoints;  % number of data points for classification

if isfield (opts, 'geneNames')
    geneNames = opts.geneNames;
else
    geneNames = string(cellstr(transpose(char( 64+(1:numRegulators ) ))));  % A-Z
end

grn.Tgg     = NaN (numGenes, numRegulators); % allocate
grn.hg      = NaN (numGenes, 1);     % allocate
grn.Rg      = NaN (numGenes, 1); % foolproofint
grn.lambdag = NaN (numGenes, 1);
grn.Dg      = NaN (numGenes, 1);

yntg = NaN (numNuclei, numTimepoints, numGenes+numExternals);
vntg = NaN (numNuclei, numTimepoints, numGenes+numExternals);


for g = 1:numGenes
    
    %======== Build xk (array of positions in HKGNBCT space)
    xkg  = reshape (xntg, numDatapoints, numGenes+numExternals);
    
    %======== Build yk (array of on/off classifications in HKGNBCT space)
    ynt = [];  % These are "group" variables (1, 0, or -1)
    vnt = []; % Storing the slope fo %estimating R, lambda, and D
    
    for n = 1:numNuclei
        mySpline = csaps (tt, xntg(n,:,g), opts.pvxOpts_ngo(n,g,1)); % Fit spline to trajectory
        mySplDer = fnder (mySpline);         % Calculate derivative w.r.t. gene g!
        for t = 1:numTimepoints              % At each timepoint
            vnt(n,t) = fnval (mySplDer, tt(t));  % Evaluate derivative
            if (abs(vnt(n,t)) > opts.pvxOpts_ngo(n,g,2))	 % If derivative is significant
                ynt(n,t) = sign (vnt(n,t));      % Infer on/off state from sign of derivative
            else							% Else infer on/off state from size of x itself
                ynt(n,t) = sign (xntg(n,t,g) - opts.pvxOpts_ngo(n,g,3));
            end           					% YLL 2018-10-2
        end
    end
    yk = reshape (ynt, numDatapoints, 1);
    
    
    yntg (:,:,g) = ynt;  % SET OPTIONAL RETURN VALUES
    vntg (:,:,g) = vnt;  % SET OPTIONAL RETURN VALUES
    
    %======== Do logistic regression (computeRegs)
    if size(unique(yk)) == 1
        warning ('Gene has only one group');
        Tg = NaN (numGenes+numExternals, 1);
        h = NaN ();
    else
        %%======== Select logistic regression implementation.
        if(strcmp(opts.lm,'glmfit'))
            % leave +1 alone, but change -1 to 0 in preparation for LogReg
            Beta = glmfit(xkg, max(yk,0), 'binomial', 'link', 'logit');  % We get h,T1,T2,...,
            Tg = Beta(2:end);
            h = Beta(1);
        elseif(strcmp(opts.lm,'FIGRlogReg'))
            Beta = FIGRlogReg(xkg, max(yk,0), opts.lambda); % We get h,T1,T2,...,
            Tg = Beta(2:end);
            h = Beta(1); 
        elseif(strcmp(opts.lm,'lassoglm'))
            [m, n] = size(xkg);
            X = xkg;
            X = [ones(m,1) X]; %% Add bias term
            B = lassoglm(X, max(yk,0)); % We get h,T1,T2,...,
            Beta = B(:,75);
            Tg = Beta(2:end);
            h = Beta(1);
        end
    end
    
    %======== Record result in Tgg array
    grn.Tgg(g,:) = Tg;
    grn.hg(g) = h;
    
    if opts.debug > 0
        fprintf(1, ['Determining R, lambda, and D for %s ' ...
            'using the %s method \n'], ...
            geneNames(g), opts.Rld_method);
    end
    
    switch opts.Rld_method
        case 'none'
            % Warning: Not inferring R, lambda, D.
            Rld_inferred = [];
        case 'slope'
            Rld_inferred(g,:) = infer_Rld_from_slope(g);
        case 'slope_nodiff'
            Rld_inferred(g,:) = infer_Rl_from_slope(g);
        case 'conc'
            Rld_inferred(g,:) = infer_Rl_from_conc(g);
        case 'kink'
            Rld_inferred(g,:) = infer_Rld_kink_approx(g);
        otherwise
            error('Unknown method for inferring R, lambda, and D %s',...
                opts.Rld_method);
    end
    
%     if opts.debug > 0
%         fprintf(1, 'Press any key to continue\n');
%         pause;
%     end  
end  % end loop over gene (g)


%======== Populate mandatory return value (struct of gene network parameters)
grn.Rg      = Rld_inferred(:,1);
grn.lambdag = Rld_inferred(:,2);
grn.Dg      = Rld_inferred(:,3);
%======== Populate optional return value (struct of diagnostic parameters)
diagnostics.yntg = yntg;
diagnostics.vntg = vntg;
return;




    function kinetic_params = infer_Rl_from_slope(g)
        
        % YLL 2019-2-4: made this simplified version
        tsafety = 3;            % this could be an option in 'opts'
        
        nmax = numNuclei;
        tmax = numTimepoints;
        gmax = numGenes+numExternals;
        yntgFilt = yntg;
        for n=1:nmax
            for t=1+tsafety:tmax-tsafety
                % FILTER OUT DATAPOINTS WHOSE SLOPES ARE
                % NOT PART OF A RUN OF 3 OR MORE OF THE SAME
                for t2=t-tsafety : t+tsafety
                    if (yntg(n,t2,g) ~= yntg(n,t,g))
                        yntgFilt(n,t,g) = 0;
                    end
                end
            end % for t
        end % for n
        
        xj = reshape (  xntg(:,:,g)'  ,[],1);  % xj = flatten (xtn)
        yj = reshape (  yntg(:,:,g)'  ,[],1);  % yj = flatten (ytn)
        yjFilt = reshape (  yntgFilt(:,:,g)'  ,[],1);
        vj = reshape (  vnt'  ,[],1);          % vj = flatten (vtn)
        
        %======== Construct the matrix Mji = [thetaj, -xj] for lsq
        thetaj = max (yj, 0);
        Mji = [thetaj, -xj];  % j=equation index, i=parameter index
        %======== Filter out datapoints close to switching events
        MjiFilt = Mji (yjFilt ~= 0, :);
        vjFilt = vj (yjFilt ~= 0);
        %======== Perform regression
        %======== This gives 3x1 col vec [R ; lambda]
        kinetic_params = lsqnonneg (MjiFilt, vjFilt);
        R      = kinetic_params(1);
        lambda = kinetic_params(2);
        D      = 0;  % populate with zeroes
        kinetic_params = [R ; lambda ; D];
        
        if opts.debug > 1
            figure(1);
            set (gcf, 'Units', 'pixels', 'Toolbar','None','MenuBar','None');
            if (g==1); clf;  end;
            
            % ASSUMES THAT grn.Tgg and grn.hg have already been populated
            % by infer()
            
            %============ COLUMN 2 ROW g: (x1, x2) SPACE
            %======== Plot concentration vs time.  YLL 2019-1-25
            subaxis (gmax,5, 2,g); hold on;
            plot (reshape(1:tmax*nmax, tmax, nmax), reshape(xj, tmax, nmax), 'k-');
            %now add the symbols            
            h1=plot (xj, 'kx');
            h2=plot (find(yjFilt>0), xj(yjFilt>0), '*', 'color', [0 .5 0]);
            h3=plot (find(yjFilt<0), xj(yjFilt<0), 'o', 'color', [1 0 0]);
            ylabel ('$x$', 'Interpreter','latex', 'Rotation', 0);
            xlim([0 tmax*nmax]);
            axis fill;
           % legend([h2 h3], 'ON', 'OFF');
            %plot trajs first with lines
            if (g==1); title('$x_g$ versus $t + t_{max} n$', 'interpreter', 'latex'); end;
            
            %============ COLUMN 3 ROW g: 
            %======== Plot "velocity" vs concentration.
            %======== Expect v = R - lambda x   for "on" points.
            %======== Expect v =   - lambda x   for "off" points.
            subaxis (gmax,5, 3,g); hold on;
            h1=plot(xj(yj > 0), vj(yj > 0), 'kx');
            h2=plot(xj(yjFilt > 0), vj(yjFilt > 0), '*', 'color', [0 .5 0]);
            h3=plot(xj(yj < 0), vj(yj < 0), 'kx');
            h4=plot(xj(yjFilt < 0), vj(yjFilt < 0), 'o', 'color', [1 0 0]);
            fplot ( @(x) (R - lambda*x), [0 2], '--', 'color', [0 .5 0]);
            fplot ( @(x) (0 - lambda*x), [0 2], '--', 'color', [1 0 0]);
            ylabel('$v$', 'interpreter', 'latex', 'Rotation', 0); 
            xlabel('$x$', 'interpreter', 'latex');
            ylim ([-2 2]);
            axis fill;
           %
           % legend([h2 h4], 'ON', 'OFF');
            if (g==1); title('$v_g$ versus $x_g$', 'interpreter', 'latex'); end;
        end
        
        
    end % from infer_Rl_from_slope()


    function kinetic_params = infer_Rld_kink_approx(g)
        
        %======== INFER the kinetic parameters (R, lambda, and D)
        %======== by modeling borders as the asymptotic (in time) solutions
        %======== of the diffusion equation with degradation. This assumes that
        %======== the gene is "on" in a domain and acts as a source for the
        %======== protein which diffuses out from the "on" domain to establish
        %======== a border.
        
        %======== We identify the timepoint with the highest expression by
        %======== smoothing with a spline and computing maximum expression.
        %======== The "on" domains are identified as intervals where
        %======== expression is greater than half max,
        %======== and the domain containing the maximum expression is chosen.
        %======== Then, for that time point, borders are identified as
        %======== consecutive points until the concentration drops below
        %======== f*Max, where Max is the maximum expression. R, l, and D are
        %======== determined by fitting the equation
        %========
        %======== x_ij = 0.5 * (R/l) * exp(-1.0*sqrt(l/D)*d_ij)
        %========
        %======== where x_ij is the concentration in the jth nucleus of the
        %======== border and d_ij is the distance from the posterior or the
        %======== anterior margin of the domain for posterior or anterior
        %======== borders respectively. For posterior borders,
        %======== d_ij = s_ij - sp_i > 0, where s_ij is the
        %======== position of the jth nucleus and sp_i is the position of the
        %======== posterior margin of the domain. For anterior borders,
        %======== d_ij = s_ij - sa_i < 0, where sa_i is the position of the
        %======== anterior margin of the domain. i denotes the domain.
        
        %======== nuclei are consecutive
        nuclei = [1:1:numNuclei]';
        
        % Loop through time points to determine when gene expression is
        % maximum.
        maxconc = 0;
        maxtimepoint = -1;
        for k=1:numTimepoints
            
            xn = xntg(:, k, g);
            spline_by_x = csaps (nuclei, xn, opts.spatialsmoothing);
            maxconc_k = -1.0*fnmin(fncmb(spline_by_x, -1.0));
            
            if (maxconc_k > maxconc)
                
                maxconc = maxconc_k;
                maxtimepoint = k;
                
            end
            
            
        end % from loop over time points
        
        %======== Identify "on" domain corresponding to the maximum expression.
        %======== A run is consecutive nuclei where concentration exceeds
        %======== half the global max. For each
        %======== such run, check whether the maximum expression within the
        %======== domain matches global maximum. If so, record the run and
        %======== stop. Otherwise, continue looking.
        xn = xntg(:, maxtimepoint, g);
        yn = ynt(:, maxtimepoint);
        
        spline_by_x = csaps (nuclei, xn, opts.spatialsmoothing);
        maxconc = -1.0*fnmin(fncmb(spline_by_x, -1.0));
        
        run_start= -1;
        run_end = -1;
        
        for j=1:numNuclei
            
            % If the point is more than half-max
            if (fnval(spline_by_x,nuclei(j)) > 0.5*maxconc)
                
                % start the run if we're not in a run
                if run_start == -1
                    run_start = j;
                end
                
                % if we're in a run and this is the last nucleus, end
                % the run
                if (run_start ~= -1) & (j == numNuclei)
                    run_end = j;
                end
                
            end
            
            % If the point is less than half-max
            if (fnval(spline_by_x,nuclei(j)) < 0.5*maxconc)
                
                % If we're in a run, then end it; do nothing if we're not
                % in a run
                if (run_start ~= -1)
                    run_end = j-1;
                end
                
            end
            
            % if the run has ended, check whether the domain contains the
            % global maximum
            if (run_end ~= -1) & (run_end > run_start+1)
                
                run_max = -1.0*fnmin(fncmb(spline_by_x, -1.0), ...
                    [nuclei(run_start),nuclei(run_end)]);
                
                if (run_max == maxconc)
                    
                    break;
                    
                end
                
                
            end % from if a run has ended
            
            % Reset the run indicators, if the run has ended
            if (run_end ~= -1)
                
                run_start = -1;
                run_end = -1;
                
            end
            
            
        end % from loop over nuclei
        
        %======== Next, we find the borders. For positions s > run_end, we
        %======== start at run_end+1 and continue until expression drops to
        %======== f*Max. For positions s < run_start, we start at run_start-1
        %======== and continue until expression drops below f*Max. The border
        %======== endpoints are saved for the next step.
        
        % A domain can have at most two borders
        border = repmat(-1, 2, 2);
        
        if run_start > 1
            j = run_start - 1;
        else
            j = 1;
        end
        
        % Find anterior (s < run_start) border
        while ((xn(j) > opts.minborder_expr_ratio*maxconc) && (j > 1))
            j = j - 1;
        end % from while to find posterior border
        
        if (run_start - j > 2)
            border(1,:) = [j,run_start-1];
        end
        
        if run_end < numNuclei
            j = run_end + 1;
        else
            j = numNuclei;
        end
        
        j = run_end + 1;
        
        % Find posterior (s > run_end) border
        while ((xn(j) > opts.minborder_expr_ratio*maxconc) && ...
                (j < numNuclei))
            j = j + 1;
        end % from while to find posterior border
        
        if (j - run_end > 2)
            border(2,:) = [run_end+1, j];
        end
        
        % Now we construct an Nx2 array, where the first column is the distance
        % of the border nucleus from the domain start or end, and the second
        % column is the concentration
        x_vs_d = [];
        
        % if the anterior border exists, record distance from domain start and
        % concentrations
        if border(1,1) ~= -1
            
            border_pos = [border(1,1):1:border(1,2)]';
            
            x_vs_d = [x_vs_d; ...
                [repmat(border(1,2)+1,length(border_pos),1) - border_pos, ...
                xn(border(1,1):border(1,2))] ];
            
        end
        
        % if the posterior border exists, record distance from domain end and
        % concentrations
        if border(2,1) ~= -1
            
            border_pos = [border(2,1):1:border(2,2)]';
            
            x_vs_d = [x_vs_d; ...
                [border_pos - repmat(border(2,1)-1,length(border_pos),1), ...
                xn(border(2,1):border(2,2))] ];
            
        end
        
        % setup the residual function, get handle to pass to lsqnonlin()
        resfunc = setResiduals(x_vs_d);
        
        [kinetic_params, resnorm] = lsqnonlin(resfunc, [10., 0.05, 0.1], ...
            [0., 0.0116, 0.], [20., 0.6931, 0.5]);
        
        if (opts.debug > 1)
            
            nucfine = [nuclei(1):0.1:nuclei(numNuclei)];
            domain = [nuclei(run_start):1:nuclei(run_end)];
            domainconcs = [repmat(0,1,run_start-1), ...
                repmat(maxconc,1,run_end-run_start+1), ...
                repmat(0,1,numNuclei-run_end)]';
            
            
            figure(6); clf; hold on;
            
            % plot data
            plot(nuclei, xn, 'k o');
            
            % plot spline fit
            plot(nucfine, fnval(spline_by_x, nucfine), 'b--');
            
            % plot "on" domains
            plot(nuclei, domainconcs, 'r--');
            
            % plot borders
            if border(1,1) ~= -1
                plot(nuclei(border(1,1):border(1,2)), ...
                    xn(border(1,1):border(1,2)), 'g*');
                
                plot(nuclei, ...
                    computeBorder(kinetic_params, nuclei(run_start), ...
                    -1, nuclei), ...
                    'k');
            end
            
            if border(2,1) ~= 1
                plot(nuclei(border(2,1):border(2,2)), ...
                    xn(border(2,1):border(2,2)), 'g*');
                
                plot(nuclei, ...
                    computeBorder(kinetic_params, nuclei(run_end), ...
                    1, nuclei), ...
                    'k');
            end
            
            ylim([0, max(max(xntg(:,:,g)))*1.1]);
            
            xlabel('AP Position');
            ylabel('x');
            
            title(sprintf('R/lambda/D inference for Gene: %s, Timepoint: %d',...
                geneNames{g}, maxtimepoint));
            
            legend({'Data', 'Spline', 'On domain', 'Border/kink tail', 'Fit'});
            
            hold off;
            
        end
        
        function fh = setResiduals(borders)
            
            %======== Returns a handle to a nested function that computes the
            %======== residuals
            %========
            %========   x_j - 0.5*(R/l)*exp(-sqrt(l/D)*d_j)
            %========
            %======== given the concentrations along borders (x_j) and the
            %======== distances (d_j) of the nuclei from the anterior or
            %======== posterior margins of the "on" domain. d_j > 0 for
            %======== posterior borders amd j_<0 for anterior borders.
            
            fh = @getResiduals;
            
            nrows = size(borders,1);
            
            % Nested function that calculates the residuals given R, l.
            function res = getResiduals(x)
                
                R = x(1);
                l = x(2);
                D = x(3);
                res = borders(:,2) - ...
                    0.5*(R/l)*exp(-1.0*sqrt(l/D)*borders(:,1));
                
            end % from getResiduals()
            
        end % from setResiduals()
        
        function concs = computeBorder(Rld, s0, border_polarity, positions)
            
            %======== Computes the "kink tail", the asymptotic approximation
            %======== to a border, given R, l, and D, the position of the kink
            %======== (s0), and the nuclei
            %========
            %======== If border is posterior, border_polarity = 1, otherwise
            %======== -1.
            
            R = Rld(1);
            l = Rld(2);
            D = Rld(3);
            
            N = length(positions);
            concs = NaN(N, 1);
            
            for j=1:N
                
                if ((positions(j) <= s0) && (border_polarity > 0)) || ...
                        ((positions(j) > s0) && (border_polarity < 0))
                    
                    
                    concs(j) = (R/l) * (1.0 - 0.5 * ...
                        exp(border_polarity*sqrt(l/D)*(positions(j)-s0)));
                    
                elseif ((positions(j) <= s0) && (border_polarity < 0)) || ...
                        ((positions(j) > s0) && (border_polarity > 0))
                    
                    concs(j) = 0.5 * (R/l) * ...
                        exp(-1.0*border_polarity*sqrt(l/D)*(positions(j)-s0));
                    
                end
                
            end % from for over positions
            
        end % from computeBorder()
        
    end % from infer_Rlk_kink_approx()

    function kinetic_params = infer_Rl_from_conc(g)
        
        %======== INFER the kinetic parameters (R and lambda), assuming that
        %======== diffusion is 0.
        
        %======== WE will identify time intervals where the gene is on.
        %======== Assuming that the rate of change when the gene is on is
        %======== R - lambda*x, where x is the concentration of the gene
        %======== product, x(t) = x(0)*exp(-l*t) + (R/l)*(1-exp(-l*t)). Let
        %======== a there be a continuous run, i, of "on" points, x_{i0},
        %======== ...,x_{ip} at times t_{i1},...,t_{ip}. This leads to p
        %======== equations of the form
        %========
        %========   x_{ij} = x_{i0}exp(-l*t_{ij}) + (R/l)*(1-exp(-l*t_{ij})),
        %========
        %======== where j=1,...,p. In an embryo there can be many such runs,
        %======== occuring over different nuclei or multiple times within
        %======== one nucleus. Let the number of such runs be r. Therefore
        %======== i=1,...,r and the total number of equations is r*p.
        %========
        %======== The algorithm, first identifies such runs of length 2 or
        %======== more, and the sets up a nonlinear optimization problem to
        %======== estimate R and lambda.
        
        % runs are stored in rows as pairs of the concentration at the starting
        % point and concentrations at subsequent points in the run.
        onruns = [];
        
        for j = 1:numNuclei
            
            run_start= -1;
            run_end = -1;
            
            for k = 1:numTimepoints
                
                % If the point is "on"
                if (ynt(j,k) > 0)
                    
                    % start the run if we're not in a run
                    if run_start == -1
                        run_start = k;
                    end
                    
                    % if we're in a run and this is the last time point, end
                    % the run
                    if (run_start ~= -1) & (k == numTimepoints)
                        run_end = k;
                    end
                    
                end
                
                % If the point is "off"
                if (ynt(j,k) < 0)
                    
                    % If we're in a run, then end it; do nothing if we're not
                    % in a run
                    if (run_start ~= -1)
                        run_end = k-1;
                    end
                    
                end
                
                % if the run has ended and has more than three points,
                % update the onruns matrix
                if (run_end ~= -1) & (run_end > run_start+1)
                    runtimepoints = tt(run_start+1:run_end) ...
                        - repmat(tt(run_start), run_end-run_start, 1);
                    runstartconc = repmat(xntg(j,run_start,g), ...
                        run_end-run_start,1);
                    runconcs = xntg(j,(run_start+1):run_end,g)';
                    
                    onruns = [onruns; [runtimepoints, runstartconc, ...
                        runconcs]];
                    
                    % Debugging plots showing the run superimposed on the "on"
                    % and "off" points
                    if (opts.debug > 1)
                        
                        drtp = runtimepoints + ...
                            repmat(tt(run_start), run_end-run_start, 1)
                        
                        figure(6); clf; hold on;
                        plot(tt(ynt(j,:) > 0), xntg(j, ynt(j, :) > 0,g), 'g*');
                        plot(tt(ynt(j,:) < 0), xntg(j, ynt(j, :) < 0,g), 'r+');
                        
                        plot(drtp, runconcs, 'bo');
                        
                        resfunc = setResiduals([runtimepoints, ...
                            runstartconc, ...
                            runconcs]);
                        [Rl, resnorm, residuals, exitflag, output] = ...
                            lsqnonlin(resfunc, ...
                            [10., 0.05], [0., 0.], [100., 1.])
                        
                        tfine = [tt(run_start):(50.-tt(run_start))/100.:50.]';
                        
                        pred_concs = computeOnTraj(Rl, runstartconc(1), ...
                            tfine-repmat(tt(run_start),101,1));
                        
                        plot(tfine, pred_concs, 'k--');
                        
                        ylim([0,250]);
                        xlim([0,60]);
                        
                        title(sprintf('Gene: %s, Nucleus: %d', ...
                            geneNames{g}, j));
                        
                        hold off;
                        
                        pause;
                        
                    end
                    
                end % from if a run has ended
                
                % Reset the run indicators, if the run has ended
                if (run_end ~= -1)
                    
                    run_start = -1;
                    run_end = -1;
                    
                end
                
                
            end % loop over time
            
        end % loop over nuclei
        
        % setup the residual function, get handle to pass to lsqnonlin()
        resfunc = setResiduals(onruns);
        
        [Rl, resnorm] = lsqnonlin(resfunc, [10., 0.05], [0., 0.], [100., 1.]);
        
        kinetic_params = [Rl, 0];
        
        function fh = setResiduals(runs)
            
            %======== Returns a handle to a nested function that computes the
            %======== residuals
            %========
            %========   x_{ij} - x_{i0}exp(-l*t_{ij})
            %========                       - (R/l)*(1-exp(-l*t_{ij})),
            %========
            %======== given the vector of concentrations (x_{ij}) observed
            %======== during a run (a contiguous set of "on" states), the
            %======== timepoints (t_{ij}), and the starting concentrations of
            %======== each run (x_{i1}). i indexes the run, and j the
            %======== timepoints within it. In the "runs" argument, the first
            %======== column contains t_{ij}, the second contains x_{i0},
            %======== and the third contains x_{ij}.
            
            fh = @getResiduals;
            
            nrows = size(runs,1);
            
            % Nested function that calculates the residuals given R, l.
            function res = getResiduals(x)
                
                R = x(1);
                l = x(2);
                exp_term = exp(-1.0*l*runs(:,1));
                R_by_l_vec = repmat(R/l, nrows, 1);
                
                
                res = runs(:,3) - ...
                    (runs(:,2) - R_by_l_vec).*exp_term - ...
                    R_by_l_vec;
                
            end % from getResiduals()
            
        end % from setResiduals()
        
        function concs = computeOnTraj(Rl, x0, times)
            
            %======== Computes the trajectory at timepoints "times" given R, l,
            %======== and the initial concentration x0, assumed to occur at t=0.
            
            R = Rl(1);
            l = Rl(2);
            
            N = length(times);
            R_by_l_vec = repmat(R/l, N, 1);
            x0_vec = repmat(x0, N, 1);
            exp_term = exp(-1.0*l*times);
            concs = (x0_vec - R_by_l_vec).*exp_term + R_by_l_vec;
            
        end % from computeOnTraj()
        
    end % from infer_Rl_from_conc()

    function kinetic_params = infer_Rld_from_slope(g)
        
        %======== INFER the kinetic parameters (R, lambda, D)
        %======== WE will solve the following set of over-determined
        %======== linear equations (perform linear least-squares regression):
        %========       s(j,tl) = R*theta(j,tl) - lambda*x(j,tl) +
        %========                   D*(x(j+1,tl)+x(j-1,tl)-2*x(j,tl)),
        %======== for j=1,...,M-1,
        %========       s(0,tl) = R*theta(0,tl) - lambda*x(0,tl) +
        %========                   D*(x(1,tl) - x(0,tl)),
        %======== for j=0, and
        %========       s(M,tl) = R*theta(M,tl) - lambda*x(M,tl) +
        %========                   D*(x(M-1,tl) - x(M,tl)).
        %======== Here, j is nucleus number, tl is the time, M is the number of
        %======== nuclei, s is the slope (not sign, but actual estimated
        %======== slope), theta is 1 if the gene is "on" and 0 if the gene is
        %======== "off", and x is the concentration of gene g.
        
        %======== Spline approximations are especially inaccurate near
        %======== points where the state switches.
        % Manu's suggested alternative: identify local
        %======== maxima and minima, and leave those and neighboring points
        %======== out.   YLL is doing something quite similar.
                
        % YLL 2019-1-30: changed algorithm and renamed many variables.
        % e.g., "s" is now "vnt" ("velocity")
        %
        % New algorithm: filter out datapoints whose "velocity" changes
        % sign between time t-tsafety and t+tsafety

        tsafety = opts.Rld_tsafety;
        
        nmax = numNuclei;
        tmax = numTimepoints;
        gmax = numGenes+numExternals;
        yntgFilt = yntg;
        for n=1:nmax
            for t=1+tsafety:tmax-tsafety
                % FILTER OUT DATAPOINTS WHOSE SLOPES ARE
                % NOT PART OF A RUN OF 3 OR MORE OF THE SAME
                for t2=t-tsafety : t+tsafety
                    if (yntg(n,t2,g) ~= yntg(n,t,g))
                        yntgFilt(n,t,g) = 0;
                    end
                end
            end % for t
        end % for n
        

        xj = reshape (  xntg(:,:,g)'  ,[],1);  % xj = flatten (xtn)
        yj = reshape (  yntg(:,:,g)'  ,[],1);  % yj = flatten (ytn)
        yjFilt = reshape (  yntgFilt(:,:,g)'  ,[],1); 
        vj = reshape (  vnt'  ,[],1);          % vj = flatten (vtn)

        %======== Construct the matrix Mji = [thetaj, -xj, ddxj]
        %======== for least squares regression
        
        %======== theta vector is just yk (after it has been mapped to (0,1)),
        %======== but we're redoing it here in case the logistic regression
        %======== gets its own function in the future.
        thetaj = max (yj, 0);
        
        %======== the weights for D can be constructed from the xntg matrix.
        %======== "dd" stands for the 2nd derivative w.r.t. space (Laplacian)
        %======== in Fick's eqn.
        ddxnt = nan(numNuclei, numTimepoints);
        
        %======== Nuclei excluding edges
        ddxnt(2:(numNuclei-1),:) = xntg(3:numNuclei,:,g) ...
            + xntg(1:(numNuclei-2),:,g) ...
            - 2.0*xntg(2:(numNuclei-1),:,g);
        
        %======== Left edge of the modeled region
        ddxnt(1,:) = xntg(2,:,g) - xntg(1,:,g);
        
        %======== Right edge of the modeled region
        ddxnt(numNuclei,:) = xntg(numNuclei-1,:,g) - xntg(numNuclei,:,g);
        
        %======== Flatten the 2D matrix into 1D vector columnwise (all nuclei
        %======== for 1st timepoint, then all nuclei for 2nd timepoint etc.)
        % YL: FLIPPED ORDERING
        ddxj = reshape (ddxnt', numDatapoints, 1);  
        
        %======== The matrix for least-squares regression
        Mji = [thetaj, -xj, ddxj];  % j=equation index, i=parameter index
        
        %======== Filter the matrix and s for noisy points
        MjiFilt = Mji (yjFilt ~= 0, :);
        vjFilt = vj (yjFilt ~= 0);
        %fprintf(1, 'Number of equations after filtering: %d\n', ...
        %    size(MjiFilt, 1));
        
        %======== Perform regression
        %======== This gives 3x1 col vec [R ; lambda ; D]
        kinetic_params = lsqnonneg (MjiFilt, vjFilt);
        R      = kinetic_params(1);
        lambda = kinetic_params(2);
        D      = kinetic_params(3);
                
        if opts.debug > 1
            figure(21);
            set (gcf, 'Units', 'pixels', 'Toolbar','None','MenuBar','None');
            if (g==1); clf;  end;
            
            %======== Plot concentration vs time.  YLL 2019-1-25
            subplot(gmax,2, g+gmax*0); hold on;            
            plot (xj, 'k+');
            plot (find(yjFilt>0), xj(yjFilt>0), 'o', 'color', [0 .5 0]);
            plot (find(yjFilt<0), xj(yjFilt<0), 'o', 'color', [1 0 0]);            
            ylabel('x'); xlabel('t + tmax*n');
            legend('All', 'Good ON', 'Good OFF');
            title (sprintf('Gene %d: x',g));
            
            %======== Plot "velocity" vs concentration.  
            %======== Expect v = R - lambda x   for "on" points.
            %======== Expect v =   - lambda x   for "off" points.
            subplot(gmax,2, g+gmax*1); hold on;
            plot(xj(yj > 0), vj(yj > 0), 'k+');
            plot(xj(yjFilt > 0), vj(yjFilt > 0), 'o', 'color', [0 .5 0]);
            plot(xj(yj < 0), vj(yj < 0), 'k+');
            plot(xj(yjFilt < 0), vj(yjFilt < 0), 'o', 'color', [1 0 0]);
            fplot ( @(x) (R - lambda*x), [0 2], '--', 'color', [0 .5 0]);
            fplot ( @(x) (0 - lambda*x), [0 2], '--', 'color', [1 0 0]);
            ylabel('v=dx/dt'); xlabel('x');
            legend('Unfiltered', 'Filtered', 'Fit', 'Fit');
            title('Slope v versus conc x');
        end
        
        
    end % from infer_Rld_from_slope()

end % from infer()

