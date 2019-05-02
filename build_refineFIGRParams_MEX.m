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
                        
                        
% generate code                
codegen -report refineFIGRParams.m
