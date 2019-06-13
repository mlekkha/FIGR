function [] = mainToyModels (runParams, debug)   % basically David's dataGen.m

% INPUTS:
%   runParams   struct with path, filenames, and unfold command line
%               options. If not provided, then will assume original paths
%               as set by Yen Lee
%   RETURNS     Nothing

%======== DEFINE default values of parameters
defaultParams.circdir = '~/genecircuits/';  %CURRENTLY NOT BEING USED PROPERLY
defaultParams.datfile = '~/genecircuits/Xntg.dat';
defaultParams.unfold = '/usr/local/fly/fly/unfold';  %HARDWIRED THE FOLDER NAME
defaultParams.timepoints = '~/genecircuits/timepoints.dat';
defaultParams.tempcirc = '~/genecircuits/networkPars';
defaultParams.unfoldout = '~/genecircuits/trajectories.dat';
defaultParams.circuit = '~/genecircuits/hkgn58c14k1_002';
defaultParams.slopethresh = 0.01;  % FOR TOY MODEL WITH TIMESCALES ABOUT 1 UNIT
defaultParams.exprthresh = 0.5;    % FOR TOY MODEL WITH 0<x<1 TYPICALLY
defaultParams.splinesmoothing = 0.99;  % FOR TOY MODEL WIT
defaultParams.unfoldopts = '-s r4 -i 0.01 -a 0.001 -x input -g h ';
defaultParams.Rld_method = 'none';        % Rld stuff not used here    --- YET!!! defaultParams.Rld_thresh_on_frac = 0.2;
defaultParams.Rld_thresh_off_frac = 0.01;
%======= SET defaults based on number of arguments provided
%======== IF no parameter struct was provided, use default struct
if nargin < 1
    runParams = defaultParams;
end
if nargin < 2
    debug = 0;
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
% rng (12345);  		% seed the random number generator reproducibly
rng ();  				% seed the random number generator irreproducibly
numGenes = 2;        	% hardwired for now
numExternals = 0;     	%

numTMatrices = 1000;	% number of different circuit parameters to simulate
numNuclei = 1000; 		% number of "nuclei" or "initial conds"

tt = (0: 0.05: 2)';
numTimepoints = numel (tt);
nuclei  = (9000+1:9000+numNuclei)';   % indices of nuclei
geneNames = {'A','B','C','D','E','F','G','H'};  % hardwired for now; could make auto gene names

TTheorList    = [];
TInferList = [];


resultsDir = sprintf ('~/genecircuits/%dgenes_%dtmats_%dnuclei/', numGenes,numTMatrices,numNuclei);
mkdir (resultsDir);
delete ([resultsDir '/TTheor.dat']); % prepare for this function to record
delete ([resultsDir '/TInfer.dat']);
delete ([resultsDir '/diagnostics.dat']);
%delete ([resultsDir '/discrep.dat']);



for m=1:numTMatrices
	fprintf ('======== Working on T-matrix no. m=%d ========\n', m);

    Rg      = ones(1,numGenes);
    lambdag = ones(1,numGenes);

    %======== GENERATE RANDOM Tgg (see David's GenTMatrices.m)
	%	Choose a point r0 within the unit hypercube using rand().
	%	Choose a unit random vector T in N=(numGenes+numExternals) dimensions,
	%		by generating a vector of N normal deviates and normalizing it.
	%	Consider the hyperplane T.(r - r0) = 0, that is, T.r + h = 0
	%		where h = -T.r0.
	
	Tgg = NaN (numGenes, numGenes+numExternals+1);
    for g=1:numGenes
		r0vec    = rand ([numGenes+numExternals 1]);
		Tvec     = normr (normrnd (0, 1, [1 numGenes+numExternals]));
		h        = -Tvec * r0vec;
		Tgg(g,:) = [Tvec h];
		
		% could put in WHILE loop
		% and keep generating Tvec until we get one
		% whose g'th component(i.e., diag element of Tmat)
		% is >0  or <0
		% depending on hether we want to study autoactivation etc.2
	end

    %======== GENERATE RANDOM INITIAL CONDITIONS Xn0g (see David's GenNFlyConcs.m)
    %	For now, choose a random point r0 within the unit hypercube using rand
    Xn0g    = rand (numNuclei, numGenes+numExternals);
    Xntg    = NaN ([numNuclei,numTimepoints,numGenes+numExternals]);
    Xntg(:,1,:) = Xn0g;
    Xntg(:,:,numGenes+1:end) = 0.;
    
    %======== computeTrajsUsingUnfold
    [xntg] = ...
    computeTrajsUsingUnfold (runParams,Rg,lambdag,Tgg,Xntg,nuclei,tt,geneNames,debug);
    
    %======== infer
    [TggInfer,Rld_inferred,diagnostics] = ...
    infer (runParams,numGenes,numExternals,xntg,nuclei,tt, geneNames,debug);
        
    %======== NORMALIZE AND COMPARE
    % YLL 2019-6-13: Don't calculate discrepancies here;
    % do it when plotting, in plotDiscrepanciesByRow.m
    for g=1:numGenes       
 	    ThTheor = normr (Tgg(g,:));        	% normalize such that ||TTheor||=1
		TTheor = ThTheor (:, 1:numGenes); 	% drop the h component
	 	TTheor = normr (TTheor);  			% normalize AGAIN
		if isnan(TggInfer(g,1)) 
	        ThInfer = NaN (1,numGenes);  
			%discrep = NaN; 					% uninferrable
		else
 	        ThInfer = normr (TggInfer(g,:));  	% normalize such that ||TInfer||=1
			TInfer = ThInfer (:, 1:numGenes); 	% drop the h component
			TInfer = normr (TInfer);			% normalize AGAIN
	 		%discrep = norm (TInfer - TTheor);   % distance between unit vectors      
		end
				
		dlmwrite ([resultsDir '/TTheor.dat'], ThTheor, 'delimiter','\t', 'precision',8,'-append');
		dlmwrite ([resultsDir '/TInfer.dat'], ThInfer, 'delimiter','\t', 'precision',8,'-append');
		%dlmwrite ([resultsDir '/discrep.dat'], discrep, 'delimiter','\t', '-append');
    end
    dlmwrite ([resultsDir '/diagnostics.dat'], diagnostics, 'delimiter','\t', '-append');
end

fprintf ('Use plotDiscrepancies([]) to visualize the results. \n');
end



