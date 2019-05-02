function grnREF = refineFIGRParams(grn, xntgFLAT, tt)
% Takes a FIGR inferred grn and refines it by local search, minimizing the
% squared difference between data (xntgEXPT) and model output. tt is the
% set of time points.  Serves as a wrapper function for conversion to C/MEX
% code. 

% declare input types

%grn
assert (isstruct(grn));

assert (isa(grn.Tgg, 'double'));
assert (all(size(grn.Tgg) > [0 0]));
assert (all(size(grn.Tgg) <= [20 20]));

assert (isa(grn.hg, 'double'));
assert (all(size(grn.hg) >= [0 1]));
assert (all(size(grn.hg) <= [20 1]));

assert (isa(grn.Rg, 'double'));
assert (all(size(grn.Rg) >= [0 1]));
assert (all(size(grn.Rg) <= [20 1]));

assert (isa(grn.lambdag, 'double'));
assert (all(size(grn.lambdag) >= [0 1]));
assert (all(size(grn.lambdag) <= [20 1]));

assert (isa(grn.Dg, 'double'));
assert (all(size(grn.Dg) >= [0 1]));
assert (all(size(grn.Dg) <= [20 1]));

%xntgFLAT
assert(isa(xntgFLAT, 'double'));
assert(all(size(xntgFLAT) > [0 0]));
assert(all(size(xntgFLAT) <= [10000 20]));

%tt
assert(isa(tt, 'double'));
assert(all(size(tt) >= [1 1]));
assert(all(size(tt) <= [100 1]));

% global options
global opts;

% "declare vectors/matrices"
init_paramvec = nan (numel(grn.Tgg) + numel(grn.hg) + numel(grn.Rg) + ...
                        numel(grn.lambdag) + numel(grn.Dg), 1);

paramREF = nan (numel(grn.Tgg) + numel(grn.hg) + numel(grn.Rg) + ...
                        numel(grn.lambdag) + numel(grn.Dg), 1);

numGenes = size(grn.Tgg, 1);

xntgEXPT = reshape(xntgFLAT, 58, 9, 7);

[init_chisq, init_rms, fh_chisq, init_paramvec] = ...
                            initChiSquare(opts, grn, xntgEXPT, tt);
                            
% optimization options
optimopts = optimset('Display', 'iter', ...
                            'MaxFunEvals', 200*length(init_paramvec), ...
                            'MaxIter', 200*length(init_paramvec));
                        
[paramREF,scoreREF,exitflag,output] = ...
                            fminsearch(fh_chisq,init_paramvec,optimopts);

grnREF = unpackParams(paramREF, grn, numGenes);

fprintf(1, 'Score: %f, Exit: %f, Count: %f\n', ...
                    scoreREF, exitflag, output.funcCount);                

end % from function declaration
