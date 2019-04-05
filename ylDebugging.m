fprintf ('\n\n');
disp ('&&&&&&&&&&&&&&&&&& TOY MODEL &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp ('&  2-gene network for illlustration purposes              &');
fprintf ('\n');

%======== User-supplied pars
 rng (12346);  		% seed the random number generator reproducibly
%rng ();  				% seed the random number generator irreproducibly
numGenes = 2;        	% hardwired for now
numTMatrices = 10;   	% number of different circuit parameters to simulate
numNuclei = 10; 		% number of "nuclei" or "initial conds"
tt = (0: 0.05: 2)';  % 0.05

%======== Derived pars
numExternals = 0;
numRegs = numGenes;
numTimepoints = numel (tt);

path (path(), '~/code/glassmodels/subaxis');

opts.debug = 0;
opts.slopethresh = .01;
opts.exprthresh = .2;
opts.splinesmoothing = 1.00;   % 1 = no smoothing, 0 = extreme smoothing
opts.Rld_method = 'slope';
opts.Rld_thresh_on_frac = 0.4;  %0.2
opts.Rld_thresh_off_frac = 0.4;  % 0.01
opts.synthesisfunction = 'synthesis_heaviside';
opts.ODEAbsTol = 1e-4;
opts.ODEsolver = 'ode45';

opts.geneNames = {'A', 'B', 'C', 'D', 'E'};

%resultsDir = sprintf ('~/genecircuits/%dgenes_%dtmats_%dnuclei/', numGenes,numTMatrices,numNuclei);
resultsDir = sprintf ('~/code/glassmodels/');
mkdir (resultsDir);
fidGRNTOY = fopen ([resultsDir '/grnTOY.dat'], 'w');
fidGRNCBI = fopen ([resultsDir '/grnCBI.dat'], 'w');
fidDISCREP = fopen ([resultsDir '/discrep.dat'], 'w');
%delete ([resultsDir '/grnTOY.dat']); % prepare for this function to record
%delete ([resultsDir '/grnCBI.dat']);
%delete ([resultsDir '/discrep.dat']);





opts.debug = 2;


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
        % could put in WHILE loop
        % and keep generating Tvec until we get one
        % whose g'th component(i.e., diag element of Tmat) is >0  or <0
        % depending on hether we want to study autoactivation etc.2
    end
    
    %======== GENERATE RANDOM KINETIC PARAMETERS R,lambda,D
    %grn.Rg      = ones (numGenes, 1);          % uniformly distributed random numbers in [0,1]
    %grn.lambdag = ones (numGenes, 1);  % uniformly distributed on [0.5, 2]
    grn.Rg      = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
    grn.lambdag = rand (numGenes, 1)*1.5+0.5;  % uniformly distributed on [0.5, 2]
    %     grn.Rg      = rand (numGenes, 1)*.1+0.9;
    %     grn.lambdag = rand (numGenes, 1)*.1+0.9;
    grn.Dg      = zeros (numGenes, 1);
    
    %======== GENERATE RANDOM INITIAL CONDITIONS x WITHIN BIOLOGICALLY RELEVANT SUBSPACE
    xntg        = NaN ([numNuclei,numTimepoints,numRegs]);
    xntg(:,1,:) = rand (numNuclei, numRegs) .* grn.Rg' ./ grn.lambdag';
    xntg(:,:,numGenes+1:end) = 0.;
    
    %======== grn ---> xntg ========
    [xntg] = computeTrajs (opts, grn, xntg, tt);
    
    disp (size(xntg));
    %======== xntg ---> grnCBI
    [grnCBI, diagnostics] = infer (opts, xntg, tt, numGenes);
    yntg = diagnostics.yntg;
    
    xkg = reshape (xntg, [numNuclei*numTimepoints numRegs]);
    ykg = reshape (yntg, [numNuclei*numTimepoints numRegs]);
    
    
    
    
    
    %======================== DISPLAY SWITCHING BOUNDARIES IN FIG 1
    %======== CALCULATE GRAPHICS DIMENSIONS
    figure (1);
    set (groot,'defaultAxesFontSize',12);
    set (groot,'defaultTextFontSize',12);
    markerSize = 14;
    set (gcf, 'Units', 'pixels', 'Position',[0 0 800 400],'Toolbar','None','MenuBar','None');
    
    for g=1:2
        yk = ykg(:,g);
        
        %======================== ACTUAL SWITCHING BOUNDARY
        subaxis (2,2,  1,g , 'Margin', 0.05); cla; hold on;
        f = @(x,y)   grn.Tgg(g,1)*x + grn.Tgg(g,2)*y + grn.hg(g);
        fcontour (f, [0 1 0 1], 'LevelList', [0], 'LineColor', 'b', 'LineWidth', 1);
        yticks ([0 .5 1]);
        xlabel ('$x_1$', 'Interpreter','latex');
        ylabel ('$x_2$', 'Interpreter','latex', 'Rotation', 0);
        if (g==1) ; xticks ([]); xlabel (''); end
        axis ([0 1 0 1]);
        %======================== INFERRED SWITCHING BOUNDARY
        f = @(x,y)   grnCBI.Tgg(g,1)*x + grnCBI.Tgg(g,2)*y + grnCBI.hg(g);
        fcontour(f, [0 1 0 1], 'LevelList', [0], 'LineStyle', ':', 'LineColor', 'b', 'LineWidth', 2); %, 'DisplayName', 'Inferred SB');
        if (g==1) ; xticks ([]); end
        %======================== ON/OFF POINTS
        negs = find (yk<=0);
        poss = find (yk>0);
        scatter (xkg(negs,1), xkg(negs,2), markerSize, [1 0 0], 'o');
        scatter (xkg(poss,1), xkg(poss,2), markerSize, [0 .8 0], '*');
        scatter (xntg(:,1,1), xntg(:,1,2), 60, [1 0 1], 'p');  % PURPLE STAR MEANS I.C.
        xlabel ('$x_1$', 'Interpreter','latex');
        ylabel ('$x_2$', 'Interpreter','latex', 'Rotation', 0);
        yticks ([0 .5 1]);
        if (g==1) ; xticks ([]); xlabel (''); end
        axis ([0 1 0 1]);
    end
    
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
    %dlmwrite ([resultsDir '/grnTOY.dat'], grnPACKED, 'delimiter','\t', 'precision',8,'-append');
    %dlmwrite ([resultsDir '/grnCBI.dat'], grnCBIPACKED, 'delimiter','\t', 'precision',8,'-append');
    writeMatrix (fidGRNTOY, grnPACKED);
    writeMatrix (fidGRNCBI, grnCBIPACKED);
    
    fprintf ('grnTRUE\n'); 
    fprintf ('%14s', 'T1g', 'T2g', 'h', 'R', 'lambda', 'D');
    fprintf ('\n'); 
    writeMatrix (1, grnPACKED);
    fprintf ('grnCBI\n'); 
    fprintf ('%14s', 'T1g', 'T2g', 'h', 'R', 'lambda', 'D');
    fprintf ('\n'); 
    writeMatrix (1, grnCBIPACKED);
    
    d1 = norm (grnCBI.Tgg(:) - grn.Tgg(:));  % flattened!
    d2 = norm (grnCBI.hg - grn.hg);
    d3 = norm (grnCBI.Rg - grn.Rg);
    d4 = norm (grnCBI.lambdag - grn.lambdag);
    d5 = norm (grnCBI.Dg - grn.Dg);
    discrep = [d1 d2 d3 d4 d5]; % row vector
    %dlmwrite ([resultsDir '/discrep.dat'], discrep, 'delimiter','\t', '-append');%dlmwrite ([resultsDir '/diagnostics.dat'], diagnostics, 'delimiter','\t', '-append');
    writeMatrix (fidDISCREP, discrep);
    
    
    %======================== DISPLAY GRNs  STILL IN FIGURE 1
    subaxis (2,2,  2,1); cla; hold on;
    rowNames = strsplit(num2str(1:numGenes));
    colNames = [strsplit(num2str(1:numRegs))  'h' 'R' '\lambda' 'D'];
    numRows = numel (rowNames);
    numCols = numel (colNames);
    axis image; axis ([.5 numCols+.5 .5 numRows+.5]);
    set(gca,'xtick',[1:numCols],'xticklabel', colNames);
    set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
    title ('grn: T_{gf}, h_g, R_g, \lambda_g, D_g');
    for r=1:numRows
        for c=1:numCols
            text (c, r, sprintf ('%0.3f', grnPACKED(r, c) ), 'HorizontalAlignment','center','FontSize',14);
        end
    end
    subaxis (2,2,  2,2); cla; hold on;
    rowNames = strsplit(num2str(1:numGenes));
    colNames = [strsplit(num2str(1:numRegs))  'h' 'R' '\lambda' 'D'];
    numRows = numel (rowNames);
    numCols = numel (colNames);
    axis image; axis ([.5 numCols+.5 .5 numRows+.5]);
    set(gca,'xtick',[1:numCols],'xticklabel', colNames);
    set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
    title ('grn: T_{gf}, h_g, R_g, \lambda_g, D_g');
    for r=1:numRows
        for c=1:numCols
            text (c, r, sprintf ('%0.3f', grnCBIPACKED(r, c) ), 'HorizontalAlignment','center','FontSize',14);
        end
    end
    
    pause;
    
end

function writeMatrix (fileID, matrix)
for r=1:size(matrix,1)
    fprintf (fileID, '%14.6f', matrix(r,:));
    fprintf (fileID, '\n');
end
end