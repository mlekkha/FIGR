function [hfig] = plotGRN (opts, grn)

%======== VISUALIZE GRN PARAMETERS USING GREEN/RED FOR +/- VALUES

geneNames = opts.geneNames;

%======== NORMALIZE EACH ROW OF THE T-MATRIX
numGenes = size (grn.Tgg, 1);
for g=1:numGenes
    fac = norm(grn.Tgg(g,:));
    grn.Tgg(g,:) = grn.Tgg(g,:) / fac;
    grn.hg(g,:)  = grn.hg(g,:) / fac;
end

mat = [grn.Tgg grn.hg grn.Rg grn.lambdag grn.Dg];
numGenes = size (grn.Tgg, 1);
numRegs  = size (grn.Tgg, 2);
rowNames = opts.geneNames(1:numGenes);
colNames = [opts.geneNames(1:numRegs)  'h' 'R' '\lambda' 'D'];
numRows = numel (rowNames);
numCols = numel (colNames);

%set(gca,'FontSize',12);
%set(0,'DefaultAxesFontSize',12);

hfig = figure();
set (hfig, 'Toolbar','None','MenuBar','None');
%set (hfig, 'pos',[0 800 800 600], 'Toolbar','None','MenuBar','None');

%myColormap = .2*[5 0 0;5 1 1;5 2 2;5 3 3;5 4 4;5 5 5;4 5 4;3 5 3;2 5 2;1 5 1;0 5 0];
myColormap = .2*[5 2 2;5 2.5 2.5; 5 3 3;5 3.5 3.5;5 4 4;5 5 5;4 5 4;3 5 3;2 5 2;1 5 1;0 5 0];
colormap (myColormap);
imagesc (mat, .05*[-1 1]);
colorbar;
axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:numCols],'xticklabel', colNames);
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
title ('grn: T_{gf}, h_g, R_g, \lambda_g, D_g');
for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat(r, c) ), 'HorizontalAlignment','center','FontSize',14);
    end
end
return;