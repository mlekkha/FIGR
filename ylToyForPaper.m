fprintf ('\n\n');
disp ('&&&&&&&&&&&&&&&&&& TOY MODEL &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp ('&  2-gene network for illlustration purposes              &');
fprintf ('\n');

path (path(), '~/lib/glassmodels/subaxis');

opts.debug = 0;
opts.slopethresh = .01;
opts.exprthresh = .2;
opts.splinesmoothing = 1.00;   % 1 = no smoothing, 0 = extreme smoothing
%opts.splinesmoothing = .999;            %%%%%%%%%%%%%;   % 1 = no smoothing, 0 = extreme smoothing
opts.Rld_method = 'slope';
opts.Rld_thresh_on_frac = 0.2;
opts.Rld_thresh_off_frac = 0.01;
%opts.synthesisfunction = 'synthesis_sigmoid_sqrt';
opts.synthesisfunction = 'synthesis_heaviside';
opts.ODEAbsTol = 1e-4;
opts.ODEsolver = 'ode45';

%======== PREPARE INITIAL CONDITIONS AND NETWORK PARAMETERS ========
% rng (12345);  % seed the random number generator predictably
grn.Tgg     = [-.1  +1 ; -1.  0];
grn.hg      = [-0.45    0.35]';
grn.Rg      = [1        2]';
grn.lambdag = grn.Rg  ;  % [0.4   .2]';
grn.Dg      = [0     0 ]';

tt   = (0: 0.05: 3)';  % time points
xn0g = [.05  .00 ;  .05 .00]; % at least two initial conds (otherwise problems with infer)

numGenes = size (grn.Tgg, 1); % assume no "externals"
numNuclei = size (xn0g, 1);
numTimepoints = numel (tt);
opts.geneNames = {'A', 'B'};

xntg        = NaN (numNuclei,numTimepoints,numGenes);
xntg(:,1,:) = xn0g;

%======== grn ---> xntg ========
[xntg] = computeTrajs (opts, grn, xntg, tt);

%======== xntg ---> grnCBI
[grnCBI, diagnostics] = infer (opts, xntg, tt, numGenes);
yntg = diagnostics.yntg;

numNuclei     = size (xntg,1);
numTimepoints = size (xntg,2);
numRegulators = size (xntg,3);
xkg = reshape (xntg, [numNuclei*numTimepoints numRegulators]);
ykg = reshape (yntg, [numNuclei*numTimepoints numRegulators]);

%======== CALCULATE GRAPHICS DIMENSIONS
%set (groot,'defaultAxesFontName','Arial');
set (groot,'defaultAxesFontSize',12);
set (groot,'defaultTextFontSize',12);
figdir = '/Users/yloh/WORK/MANU/GRN_ms/';
markerSize = 18;
% IN CENTIMETERS
widthFigure = 18; % 19?
widthMargin = 1;
numPanelsX = 3;
numPanelsY = 2; % calculation is not really correct
widthPanel = widthFigure/numPanelsX; %(widthFigure - (numPanelsX-1)*widthMargin) / numPanelsX;
heightPanel = widthPanel * 1; % square panels
heightFigure = heightPanel * numPanelsY;
set (gcf, 'Units', 'centimeters', 'Position',[0 0 widthFigure 16],'Toolbar','None','MenuBar','None');

for g=1:2
    xt = xntg(1,:,g)';  % column vector
    yt = yntg(1,:,g)';  % column vector
    yk = ykg(:,g);
    
    %======================== 1: ACTUAL SWITCHING BOUNDARY
    subaxis (3,3,  1,g+1 , 'Spacing', 0.07, 'ML', 0.07, 'MR', 0.015, 'MT', 0.02); cla; hold on;
    
    f = @(x,y)   grn.Tgg(g,1)*x + grn.Tgg(g,2)*y + grn.hg(g);
    fcontour (f, [0 1 0 1], 'LevelList', [0], 'LineColor', 'b', 'LineWidth', 1, 'Fill', 'on');
    colormap ([1 .6 .6; .6 1 .6]);        %red or green areas
    if g==1
        myLeg = legend ('$T_{11} x_1 + T_{12} x_2 + h_1 = 0$', 'AutoUpdate', 'off' );
    else
        myLeg = legend ('$T_{21} x_1 + T_{22} x_2 + h_2 = 0$', 'AutoUpdate', 'off' );
    end
    set (myLeg, 'interpreter', 'latex');
    set (myLeg, 'location', 'best'); %[.1 .9 .9 .1]);
    
    plot (xntg(1,:,1), xntg(1,:,2), 'Color', [1 0 1]);  %traj
    
    yticks ([0 .5 1]);
    xlabel ('$x_1$', 'Interpreter','latex');
    yhandle = ylabel ('$x_2$', 'Interpreter','latex', 'Rotation', 0);
    set (yhandle, 'Units', 'Normalized', 'Position', [-0.18, 0.45, 0]);
    if (g==1) ; xticks ([]); xlabel (''); end
    axis ([0 1 0 1]);
    
    text(-0.1,1.1,char(65+3*g), 'FontSize', 12, 'FontWeight', 'bold');
    %======================== 2: TRAJECTORIES x_g(t)
    subaxis (3,3,  2,g+1 , 'Spacing', 0.07, 'ML', 0.07, 'MR', 0.015, 'MT', 0.02); cla; hold on;
    negs = find (yt<=0);
    poss = find (yt>0);
    scatter (tt(negs), xt(negs), markerSize, [1 0 0], 'o');
    scatter (tt(poss), xt(poss), markerSize, [0 .8 0], '*');
    ylim ([0 1]);
    myLeg = legend ( ...
        sprintf('$y_%d = -1$',g), ...
        sprintf('$y_%d = +1$',g), ...
        'AutoUpdate', 'off' );
    set (myLeg, 'interpreter', 'latex');
    set (myLeg, 'location', 'best');    
    
    tFineGrid = linspace (min(tt), max(tt), 501);
    mySpline = csaps (tt, xt, opts.splinesmoothing);   % WE CHECKED THAT WE NEED THIS MUCH SMOOTHING
    mySplDer = fnder (mySpline);         % Calculate derivative w.r.t. gene g!
    plot (tFineGrid, fnval (mySpline, tFineGrid), 'DisplayName', '' );
    yticks ([0 .5 1]);
    xlabel ('$t$', 'Interpreter','latex');
    if (g==1) ; xticks ([]); xlabel (''); end
    yhandle = ylabel (sprintf('$x_%d$',g), 'Interpreter','latex', 'Rotation', 0);
    set (yhandle, 'Units', 'Normalized', 'Position', [-0.18, 0.45, 0]);
      
    text(-0.33,1.1,char(66+3*g), 'FontSize', 12, 'FontWeight', 'bold');
    %======================== 3: INFERRED SWITCHING BOUNDARY
    subaxis (3,3,  3,g+1 , 'Spacing', 0.07, 'ML', 0.07, 'MR', 0.015, 'MT', 0.02); cla; hold on;
    
    f = @(x,y)   grnCBI.Tgg(g,1)*x + grnCBI.Tgg(g,2)*y + grnCBI.hg(g);
    fcontour(f, [0 1 0 1], 'LevelList', [0], 'LineStyle', '--', 'LineColor', 'b', 'LineWidth', 1); %, 'DisplayName', 'Inferred SB');
    if (g==1) ; xticks ([]); end
    if g==1
        myLeg = legend ('$\widetilde{T}_{11} x_1 + \widetilde{T}_{12} x_2 + \widetilde{h}_1 = 0$' );
    else
        myLeg = legend ('$\widetilde{T}_{21} x_1 + \widetilde{T}_{22} x_2 + \widetilde{h}_2 = 0$' );
    end
    set (myLeg, 'interpreter', 'latex');
    set (myLeg, 'location', 'best');        
    legend ('AutoUpdate', 'off');
    
    negs = find (yk<=0);
    poss = find (yk>0);
    scatter (xkg(negs,1), xkg(negs,2), markerSize, [1 0 0], 'o', 'DisplayName', [opts.geneNames{g} ' off']);
    scatter (xkg(poss,1), xkg(poss,2), markerSize, [0 0.8 0], '*', 'DisplayName', [opts.geneNames{g} ' on']);
    xlabel ('$x_1$', 'Interpreter','latex');
    yhandle = ylabel ('$x_2$', 'Interpreter','latex', 'Rotation', 0);
    set (yhandle, 'Units', 'Normalized', 'Position', [-0.18, 0.45, 0]);
    yticks ([0 .5 1]);
    if (g==1) ; xticks ([]); xlabel (''); end
    axis ([0 1 0 1]);
    
    text(-0.1,1.1,char(67+3*g), 'FontSize', 12, 'FontWeight', 'bold');
end

%======= Print actual and inferred GRNs

set(0,'DefaultAxesFontSize',10);
myhsv = [[ones(128,1); 0.4*ones(128,1)] 0.5*[[1:-1/127:0]'; [0:1/127:1]'] 0.9*ones(256,1)];
myColormap = hsv2rgb(myhsv);

%======== The CBI parameters after refinement
grntoplot = grn;

%======== NORMALIZE EACH ROW OF THE T-MATRIX
numGenes = size (grntoplot.Tgg, 1);
for g=1:numGenes
    fac = norm(grntoplot.Tgg(g,:));
    grntoplot.Tgg(g,:) = grntoplot.Tgg(g,:) / fac;
    grntoplot.hg(g,:)  = grntoplot.hg(g,:) / fac;
end

mat = [grntoplot.Tgg grntoplot.hg grntoplot.Rg grntoplot.lambdag grntoplot.Dg];
numGenes = size (grntoplot.Tgg, 1);
numRegs  = size (grntoplot.Tgg, 2);
rowNames = opts.geneNames(1:numGenes);
colNames = [opts.geneNames(1:numRegs)  'h' 'R' '\lambda' 'D'];
numRows = numel (rowNames);
numCols = numel (colNames);

subaxis (3,3,  1,1 , 'Spacing', 0.07, 'ML', 0.07, 'MR', 0.015, 'MT', 0.02); 
matcolors = [grntoplot.Tgg grntoplot.hg zeros(numGenes,3)];
im = imagesc (matcolors, 0.2*[-1 1]);
axis image; % force 1:1 aspect ratio
colorbar('Location', 'southoutside');
colormap (gca, myColormap);

% aborted attemt to have a color bar on the right of the image, same as
% image size
%cbh = colorbar(gca, 'Location', 'manual', ...
%       'Position', [oldpos(1)+oldpos(3)+0.01 (4/6)*oldpos(3) 0.01 (2/6)*oldpos(3)]);

set(gca,'xtick',[1:numCols],'xticklabel', colNames);
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
set(gca,'XAxisLocation', 'top');

%Reduce the height of the panel
oldpos = get(gca, 'Position');
heightdiff = (1/6)*oldpos(3); 
set(gca, 'Position', [oldpos(1) oldpos(2)-heightdiff oldpos(3) oldpos(4)]);

for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat(r, c) ), 'HorizontalAlignment','center','FontSize',8);
    end
end



set (gca, 'FontSize', 10);
set (gca, 'FontWeight', 'bold');
set (gca, 'TickLength', [0 0.025]);
text(-0.05,0.05,'A', 'FontSize', 12, 'FontWeight', 'bold');
%======== The inferred parameters 
grntoplot = grnCBI;

%======== NORMALIZE EACH ROW OF THE T-MATRIX
numGenes = size (grntoplot.Tgg, 1);
for g=1:numGenes
    fac = norm(grntoplot.Tgg(g,:));
    grntoplot.Tgg(g,:) = grntoplot.Tgg(g,:) / fac;
    grntoplot.hg(g,:)  = grntoplot.hg(g,:) / fac;
end

mat = [grntoplot.Tgg grntoplot.hg grntoplot.Rg grntoplot.lambdag grntoplot.Dg];
numGenes = size (grntoplot.Tgg, 1);
numRegs  = size (grntoplot.Tgg, 2);
rowNames = opts.geneNames(1:numGenes);
colNames = [opts.geneNames(1:numRegs)  'h' 'R' '\lambda' 'D'];
numRows = numel (rowNames);
numCols = numel (colNames);

subaxis (3,3,  3,1 , 'Spacing', 0.07, 'ML', 0.07, 'MR', 0.015, 'MT', 0.02); 
matcolors = [grntoplot.Tgg grntoplot.hg zeros(numGenes,3)];
imagesc (matcolors, .2*[-1 1]);
axis image; % force 1:1 aspect ratio
colorbar('Location', 'southoutside');
colormap (gca, myColormap);
set(gca,'xtick',[1:numCols],'xticklabel', colNames);
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
set(gca,'XAxisLocation', 'top');

%Reduce the height of the panel
oldpos = get(gca, 'Position');
heightdiff = (1/6)*oldpos(3); 
set(gca, 'Position', [oldpos(1) oldpos(2)-heightdiff oldpos(3) oldpos(4)]);

for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat(r, c) ), 'HorizontalAlignment','center','FontSize',8);
    end
end



set (gca, 'FontSize', 10);
set (gca, 'FontWeight', 'bold');
set (gca, 'TickLength', [0 0.025]);
text(-0.05,0.05,'C', 'FontSize', 12, 'FontWeight', 'bold');






export_fig -transparent toyNinePanels.pdf;









%======================= WRITE TO PDF FILE ======================
function printpdf(h,outfilename)
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);
end
