function [init_chisq, init_rms, fh_chisq, init_paramvec] = ...
                            initChiSquare(opts, grn, xntgEXPT, tt)

    % Returns handle to nested function that computes the value of the
    % objective function for a given GRN by solving the gene circuit 
    % equations
    %       
    %    d/dt x_ntg = R_g g( sum_f T_gf x_ntf + h_g) - lambda_g x_ntg
    %                       + D_g(x_(n+1)tg + x_(n-1)tg - 2 x_ntg)
    %
    % and computing the sum of squared differences
    %
    %    chi_sq = sum_{j=1}^{numTimepoints} sum_{k=1}{numNuclei}
    %             sum_{l=1}{numGenes} (xntgEXPT(j,k) - xntgMODEL(j,k))^2,
    %
    % where xntgEXPT are the experimentally observed gene expression
    % concentrations and xntgMODEL is the output of the model with
    % parameters specified in grn.
    %
    % InitChiSquare returns:
    % 
    %   init_chisq: the chisq for the values of parameters given in grn
    %   fh_chisq: function handle to getChiSquare()
    %   init_paramvec: grn parameters packed into a parameter vector, to be
    %       used as initial parameter set in optimization functions

    
    % number of genes being modeled
    numGenes = numel(grn.Rg);

    % number of time points
    numTimepoints = numel(tt);

    % number of nuclei
    numNuclei = size(xntgEXPT, 1);

    % number of regulators. A row of Tgg is (T_gf, M_g, E_ge, h_g), where
    % f=1,...,numGenes. M_g and E_ge represent the interconnection
    % strengths of external regulators (totalling numExternals). h_g is the
    % threshold.
    numRegulators = size(grn.Tgg, 2);
    numExternals = numRegulators - numGenes;

    % local variables for the nested function. 
    lopts = opts;
    lxntgEXPT = xntgEXPT;
    ltt = tt;

    % pack the parameter structure into an initial parameter vector
    init_paramvec = packParams(grn, numGenes);

    % return the function handle to getChiSquare
    fh_chisq = @getChiSquare;

    % determine the chi square at parameter values specified in grn
    init_chisq = fh_chisq(init_paramvec);
    init_rms = sqrt(init_chisq/((numTimepoints-1)*numNuclei*numGenes));


    % Nested function that computes the ChiSquare, given a set of GRN
    % parameters.

    function chisq = getChiSquare(paramvec)

        % Declare the grn structure and elements
        
        lgrn.Tgg = nan (numGenes, numRegulators);
        lgrn.hg = nan (numGenes, 1);
        lgrn.Rg = nan (numGenes, 1);
        lgrn.lambdag = nan (numGenes, 1);
        lgrn.Dg = nan (numGenes, 1);

        % unpack the parameter vector into the lgrn struct
        lgrn = unpackParams(paramvec, lgrn, numGenes);
                
        % Call computeTrajs to solve the equations with given values of
        % parameters
        [xntgMODEL] = computeTrajs (lopts, lgrn, lxntgEXPT, ltt);


        % Compute chi-square. Exclude the first timepoint, since it is the
        % initial condition
        diff = xntgEXPT(:,2:end,1:numGenes) - xntgMODEL(:,2:end,1:numGenes);
        
        chisq = sum(reshape(diff.*diff, ...
                                numGenes*numNuclei*(numTimepoints-1), 1));

    end % from getChiSquare()

end % from InitChiSquare()

    % saving old code in case the pack/unpack functions have problems

    %packing code
    %pos = 0;
    %init_paramvec(pos+1:pos+numGenes,1) = lgrn.Rg;

    %pos = pos + numGenes;
    %init_paramvec(pos+1:pos+numGenes*numRegulators,1) = ...
    %                reshape(lgrn.Tgg, numGenes*numRegulators, 1);

    %pos = pos + numGenes*numRegulators;
    %init_paramvec(pos+1:pos+numGenes,1) = lgrn.hg;

    %pos = pos + numGenes;
    %init_paramvec(pos+1:pos+numGenes,1) = lgrn.Dg;

    %pos = pos + numGenes;
    %init_paramvec(pos+1:pos+numGenes,1) = lgrn.lambdag;

    % unpacking code
    %pos = 0;
    %        lgrn.Rg = paramvec(pos+1:pos+numGenes);

    %        pos = pos + numGenes;
    %        lgrn.Tgg = reshape(paramvec(pos+1:pos+numGenes*numRegulators), ...
    %                            numGenes, numRegulators);

    %        pos = pos + numGenes*numRegulators;
    %        lgrn.hg = paramvec(pos+1:pos+numGenes);
    %
    %        pos = pos + numGenes;
    %        lgrn.Dg = paramvec(pos+1:pos+numGenes);

    %        pos = pos + numGenes;
    %        lgrn.lambdag = paramvec(pos+1:pos+numGenes);


