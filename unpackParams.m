function grn = unpackParams(paramvec, grn, numGenes)

    % Unpacks the parameter structure from a 1D parameter vector
    % Note that the structure grn must already exist. This is to avoid
    % unnecessary memory operations.

    % numRegulators is the number of columns of Tgg
    numRegulators = length(paramvec)/numGenes - 4;

    % reshape 1D parameter vector in matric with numGenes rows
    grnmatrix = reshape(paramvec, numGenes, numRegulators+4);

    % Read off the parameters as columns of the matrix
    grn.Rg = grnmatrix(:,1);
    grn.Tgg = grnmatrix(:,2:numRegulators+1);
    grn.hg = grnmatrix(:,numRegulators+2);
    grn.lambdag = grnmatrix(:,numRegulators+3);
    grn.Dg = grnmatrix(:,numRegulators+4);

end % from unpackParams()
