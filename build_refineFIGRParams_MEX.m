% Script build_refineFIGRParams_MEX
%
% Builds the MEX file for refineFIGRParams.m. Since opts and ODEopts are
% global variables, they are defined here first before compilation. If any
% fields are added, removed, or renamed in the rest of the code, they must
% be changed here as well. The values of the fields are immaterial, but the
% type, size, and their order of appearance in opts matter for code
% generation.

% declare global opts struct
global opts;
opts = struct(  'debug', 0, ...
                'slopethresh', 1, ...
                'exprthresh', 100.0, ...
                'splinesmoothing', 0.01, ...
                'spatialsmoothing', 0.5, ...
                'minborder_expr_ratio', 0.01, ...
                'Rld_method', 'kink', ...
                'Rld_tsafety', 3, ...
                'synthesisfunction', 'synthesis_sigmoid_sqrt', ...
                'ODEAbsTol', 1e-3, ...
                'ODEsolver', 'ode45');
            
% declare global ODEopts struct
global ODEopts;
ODEopts = odeset('AbsTol', 1e-3);
               
% declare optimization options for refinement
global optimopts;
optimopts = optimset( 'Display', 'Iter', ...
                      'MaxFunEvals', 20000, ...
                      'MaxIter', 20000);
                  
% JEH 2021-1-5
% Defining upper-boudns for the elements of the grnFIGR struct.
numGenesMax = 1000;
numNucleiMax = 1000;
numTimepointsMax = 1000;
numExternalGenesMax = 1000;

% This command uses coder.typeof to specify struct as input for 
% function refineFIGRParams. The first argument, 0, indicates the input data type (double) 
% and complexity (real). The second argument, [num* 1], indicates the size, a 
% matrix with two dimensions. The third argument, 1, indicates that the input 
% is variable sized. The upper bound is num* for the first dimension and 1 for 
% the second dimension. 
grnFIGRarg = coder.typeof( ...
struct('Tgg', coder.typeof(0, [numGenesMax numGenesMax], [1 1]), ...
       'hg', coder.typeof(0, [numGenesMax 1], [1 0]), ...
       'Rg', coder.typeof(0, [numGenesMax 1], [1 0]), ...
       'lambdag', coder.typeof(0, [numGenesMax 1], [1 0]), ...
       'Dg', coder.typeof(0, [numGenesMax 1], [1 0])));
                 
                
% JEH 2021-1-5

% User should adjust the below variables depending on the investigating 
% problem.

numGenes = 4; 
numNuclei = 58;
numTimepoints = 9;
numExternalGenes = 3;

% JEH 2021-1-6
% Wasn't able to assign upper bounds for the below variables due to this
% error: For code generation, FMINSEARCH requires FUN to return a 
% fixed-size scalar.

xntgFLATarg = coder.typeof( 0, [numNuclei*numTimepoints numGenes+numExternalGenes], 0);
 ttArg = coder.typeof(0, [numTimepoints 1], 0);

 
numGenesArg = coder.typeof(0);
numNucleiArg = coder.typeof(0);

% generate code             
codegen -report refineFIGRParams -args {grnFIGRarg, xntgFLATarg, ttArg, numGenesArg, numNucleiArg}        
