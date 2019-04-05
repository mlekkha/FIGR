function [] = runToyModel (runParams, debug)

% INPUTS:
%   runParams   struct with path, filenames, and unfold command line
%               options. If not provided, then will assume original paths
%               as set by Yen Lee
%   RETURNS     Nothing

%======== DEFINE default values of parameters
defaultParams.circdir = '~/genecircuits/';  %CURRENTLY NOT BEING USED PROPERLY
defaultParams.datfile = '~/genecircuits/Xntg.dat';
defaultParams.unfold = '/usr/local/fly/fly/unfold';  %HARDWIRED THE FOLDER NAME
defaultParams.timepoints = '~/genecircuits/timepoints';
defaultParams.tempcirc = '~/genecircuits/networkPars';
defaultParams.unfoldout = '~/genecircuits/trajectories';
defaultParams.circuit = '~/genecircuits/hkgn58c14k1_002';
defaultParams.TggTheorFile = 'TggManu.dat';
defaultParams.TggInferFile = 'TggInfer.dat';
defaultParams.slopethresh = 0.01;  % FOR TOY MODEL WITH TIMESCALES ABOUT 1 UNIT
defaultParams.exprthresh = 0.5;    % FOR TOY MODEL WITH 0<x<1 TYPICALLY
defaultParams.splinesmoothing = 0.99;  % FOR TOY MODEL WIT
defaultParams.unfoldopts = '-s r4 -i 0.01 -a 0.001 -x input -g h ';
defaultParams.Rld_method = 'conc'; %'slope';
defaultParams.Rld_thresh_on_frac = 0.2;
defaultParams.Rld_thresh_off_frac = 0.01;
%======= SET defaults based on number of arguments provided
%======== IF no parameter struct was provided, use default struct
if nargin < 1
    runParams = defaultParams;
end
if nargin < 2
    debug = 0
end
%======== PARSE the fields in runParams argument and add any missing
%======== ones with default values
param_names = fields(defaultParams);
for (j = 1:length(param_names))   
    if (~isfield(runParams, param_names{j}))        
        runParams.(param_names{j}) = defaultParams.(param_names{j});
    end
end
%======== unfold options include the timepoints file with -j
runParams.unfoldopts = [runParams.unfoldopts ' -j ' runParams.timepoints];




%======== PREPARE INITIAL CONDITIONS AND NETWORK PARAMETERS ========
% rng (12345);  % seed the random number generator predictably
numGenes = 2;   % hardwired for now
numExternals = 0;   % 1
numNuclei = 10;
tt = (0: 0.05: 2)';
numTimepoints = numel (tt);
nuclei  = (9000+1:9000+numNuclei)';    % indices of nuclei
geneNames = {'A', 'B', 'Z'};  %hardwired for now

Rg      = [1   1 ];
lambdag = [1   1 ];
Tgg     = [0   +1   -0.5;   -1   0   +0.5];   % include h_g as last column
Xn0g    = rand (numNuclei, numGenes+numExternals);
Xntg    = repmat (nan, [numNuclei,numTimepoints,numGenes+numExternals]);
Xntg(:,1,:) = Xn0g;
Xntg(:,:,numGenes+1:end) = 0.;

[xntg] = computeTrajsUsingUnfold (runParams,Rg,lambdag,Tgg,Xntg,nuclei,tt,geneNames,debug);

%======== INFER TggInfer
[TggInfer,Rld_inferred,diagnostics]  ...
= infer (runParams, numGenes,numExternals,xntg,nuclei,tt, geneNames,debug);

%======== NORMALIZE AND COMPARE
for g=1:numGenes
    Tgg(g,:) = Tgg(g,:) / norm(Tgg(g,:));
    TggInfer(g,:) = TggInfer(g,:) / norm(TggInfer(g,:));
end
Tgg = Tgg * 1000;
TggInfer = TggInfer * 1000;


  %======== Print R, l, D
    fprintf(1,'\n    \t R \t Half-life \t Max Expression \t D \n');
    for g=1:numGenes

        fprintf(1, '%s:\t %f \t %f \t %f \t %f \n', geneNames{g}, ...
                                Rld_inferred(g,1), ...
                                log(2)/Rld_inferred(g,2), ...
                                Rld_inferred(g,1)/Rld_inferred(g,2), ...
                                Rld_inferred(g,3));
    end


%======== COLOR PLOT(S)
figure (5);
plotTMatrices (Tgg, TggInfer, geneNames);






return;


%dlmwrite (runParams.TggTheorFile, Tgg, 'delimiter', '\t', 'precision', 2);
%dlmwrite (runParams.TggInferFile, TggInfer, 'delimiter', '\t', 'precision', 2);
%return;
end



