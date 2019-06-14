fprintf ('\n\n');
disp ('&&&&&&&&&&&&&&&&&& TOY MODEL &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp ('  Run inference on many 2-gene networks                             ');
disp ('  Write grnTOY.dat, grnCBI.dat, discrep.dat in appropriate directory');
disp (' ');
disp (' ');
disp ('  For 2 genes, 1000 nuclei, each inference takes about 1 minute ');
disp (' ');
disp ('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
fprintf ('\n');

%======== User-supplied pars
rng (12345);            % seed the random number generator reproducibly
%rng ();  				% seed the random number generator irreproducibly
numGenes     = 10;     	% hardwired for now
numTMatrices = 100;  	% number of different circuit parameters to simulate
numNuclei    = 30;     % number of "nuclei" or "initial conds"
tt = (0: 0.05: 2)';     % timepoints

%======== Derived pars
numExternals = 0;
numRegs = numGenes;
numTimepoints = numel (tt);

opts.debug = 0;
opts.slopethresh = .01;
opts.exprthresh = .2;
opts.splinesmoothing = 1.00;   % 1 = no smoothing, 0 = extreme smoothing
opts.Rld_method = 'slope_nodiff'; % see documentation in infer.m
opts.synthesisfunction = 'synthesis_heaviside';
opts.ODEAbsTol = 1e-4;
opts.ODEsolver = 'ode45';

resultsDir = sprintf ('%dgenes_%dtmats_%dnuclei/', numGenes,numTMatrices,numNuclei);
mkdir (resultsDir);
fidGRNTOY = fopen ([resultsDir '/grnTOY.dat'], 'w');
fidGRNCBI = fopen ([resultsDir '/grnCBI.dat'], 'w');
fidDISCREP = fopen ([resultsDir '/discrep.dat'], 'w');

for m=1:numTMatrices
    fprintf ('======== Working on T-matrix no. m=%d ========\n', m);
    
    %======== GENERATE RANDOM REGULATORY PARAMETERS T,h
    %	Choose a point r0 within the unit hypercube using rand().
    %	Choose a unit random vector T in N=(numGenes+numExternals) dimensions,
    %		by generating a vector of N normal deviates and normalizing it.
    %	Consider the hyperplane T.(r - r0) = 0, that is, T.r + h = 0
    %		where h = -T.r0.
    grn.Tgg = NaN (numGenes, numRegs);
    grn.hg = NaN (numGenes, 1);
    for g=1:numGenes
        r0vec    = rand ([numRegs 1]);
        Tvec     = normrnd (0, 1, [numRegs 1]);
        Tvec     = Tvec ./ norm(Tvec);
        h        = -Tvec' * r0vec;
        grn.Tgg(g,:) = Tvec';
        grn.hg(g)    = h;
    end   % The above could be modified so that T is auto-activating or auto-repressing
    
    %======== GENERATE RANDOM KINETIC PARAMETERS R,lambda,D
    grn.Rg      = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
    grn.lambdag = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
    grn.Dg      = zeros (numGenes, 1);
    
    %======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
    xntg        = NaN ([numNuclei,numTimepoints,numRegs]);
    xntg(:,1,:) = rand (numNuclei, numRegs) .* grn.Rg' ./ grn.lambdag';
    xntg(:,:,numGenes+1:end) = 0.;
    
    %======== grn ---> xntg ========
    [xntg] = computeTrajs (opts, grn, xntg, tt);
    
    %======== xntg ---> grnCBI
    [grnCBI, diagnostics] = infer (opts, xntg, tt, numGenes);
    
    %======== NORMALIZE AND COMPARE
    for g=1:numGenes
        mynorm          = norm( grn.Tgg(g,:) );
        grn.Tgg(g,:)    = grn.Tgg(g,:) / mynorm;
        grn.hg(g)       = grn.hg(g)    / mynorm;
        mynorm          = norm( grnCBI.Tgg(g,:) );
        grnCBI.Tgg(g,:) = grnCBI.Tgg(g,:) / mynorm;
        grnCBI.hg(g)    = grnCBI.hg(g)    / mynorm;
    end
    grnPACKED = [grn.Tgg grn.hg grn.Rg grn.lambdag grn.Dg];
    grnCBIPACKED = [grnCBI.Tgg grnCBI.hg grnCBI.Rg grnCBI.lambdag grnCBI.Dg];
    writeMatrix (fidGRNTOY, grnPACKED);
    writeMatrix (fidGRNCBI, grnCBIPACKED);
end

function writeMatrix (fileID, matrix)
for r=1:size(matrix,1)
    fprintf (fileID, '%14.6f', matrix(r,:));
    fprintf (fileID, '\n');
end
end