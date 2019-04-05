function paramvec = packParams(grn, numGenes)

    % Packs the parameter structure into a 1D parameter vector

    % numRegulators is the number of columns of Tgg
    numRegulators = size(grn.Tgg,2);

    paramvec = reshape([grn.Rg grn.Tgg grn.hg grn.lambdag ...
                            grn.Dg], numGenes*(numRegulators+4), 1);

end % from packParams()
