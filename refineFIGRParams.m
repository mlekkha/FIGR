function grn = refineFIGRParams(grn, xntgFLAT, tt, numGenes, numNuclei) 

%#codegen

% Takes a FIGR inferred grn and refines it by local search, minimizing the
% squared difference between data (xntgEXPT) and model output. tt is the
% set of time points.  Serves as a wrapper function for conversion to C/MEX
% code. 

% JEH 2021-1-6 

% Define maximum size for input parameters. 
numGenesMax = 1000;
numNucleiMax = 1000;
numTimepointsMax = 1000;
numExternalGenesMax = 1000;

% declare input types

assert (isstruct(grn));

assert (isa(grn.Tgg, 'double'));
assert (all(size(grn.Tgg) > [0 0]));
assert (all(size(grn.Tgg) <= [numGenesMax numGenesMax]));

assert (isa(grn.hg, 'double'));
assert (all(size(grn.hg) >= [0 1]));
assert (all(size(grn.hg) <= [numGenesMax 1]));

assert (isa(grn.Rg, 'double'));
assert (all(size(grn.Rg) >= [0 1]));
assert (all(size(grn.Rg) <= [numGenesMax 1]));

assert (isa(grn.lambdag, 'double'));
assert (all(size(grn.lambdag) >= [0 1]));
assert (all(size(grn.lambdag) <= [numGenesMax 1]));

assert (isa(grn.Dg, 'double'));
assert (all(size(grn.Dg) >= [0 1]));
assert (all(size(grn.Dg) <= [numGenesMax 1]));

%xntgFLAT
assert(isa(xntgFLAT, 'double'));
assert(all(size(xntgFLAT) > [0 0]));
assert(all(size(xntgFLAT) <= [numTimepointsMax numTimepointsMax]));

%tt
assert(isa(tt, 'double'));
assert(all(size(tt) >= [1 1]));
assert(all(size(tt) <= [numTimepointsMax 1]));

% global options
global opts;

% declare optimization options for refinement
global optimopts;

% "declare vectors/matrices"
init_paramvec = nan (numel(grn.Tgg) + numel(grn.hg) + numel(grn.Rg) + ...
                        numel(grn.lambdag) + numel(grn.Dg), 1);

paramREF = nan (numel(grn.Tgg) + numel(grn.hg) + numel(grn.Rg) + ...
                        numel(grn.lambdag) + numel(grn.Dg), 1);

% JEH 2021-1-5

% Instead of hard-coded values, the reshape() is performed with variables
% taken as params in refineFIGRParams().
% numGenes in reshape() below is the sum of genes and the external 
% regulators.

%xntgEXPT = reshape(xntgFLAT, 58, 9, 7);
xntgEXPT = reshape(xntgFLAT, numNuclei, size(tt,1), numGenes);

[init_chisq, init_rms, fh_chisq, init_paramvec] = ...
                            initChiSquare(opts, grn, xntgEXPT, tt);
                            
[paramREF,scoreREF,exitflag,output] = ...
                            fminsearch(fh_chisq,init_paramvec,optimopts);

numGenes = size(grn.Tgg, 1);

grnREF = unpackParams(paramREF, grn, numGenes);

fprintf(1, 'Score: %f, Exit: %f, Count: %f\n', ...
                    scoreREF, exitflag, output.funcCount);   


end % from function declaration
