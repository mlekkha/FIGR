function [] = plotGRNComparison (opts, grn1, grn2)

%======== COMPARE TWO GRNs IN A VISUALLY REVEALING WAY
% TO DO: REPLACE THIS BY A SINGLE plotGRN FUNCTION
% THAT CAN BE CALLED TWICE, AND FIGURE WINDOWS CAN BE COMPARED

geneNames = opts.geneNames;


%======== NORMALIZE EACH ROW OF THE T-MATRIX
numGenes = size (grn1.Tgg, 1);
for g=1:numGenes
    fac = norm(grn1.Tgg(g,:));
    grn1.Tgg(g,:) = grn1.Tgg(g,:) / fac;
    grn1.hg(g,:)  = grn1.hg(g,:) / fac;
    fac = norm(grn2.Tgg(g,:));
    grn2.Tgg(g,:) = grn2.Tgg(g,:) / fac;
    grn2.hg(g,:)  = grn2.hg(g,:) / fac;
end


mat1 = [grn1.Tgg grn1.hg grn1.Rg grn1.lambdag grn1.Dg];
mat2 = [grn2.Tgg grn2.hg grn2.Rg grn2.lambdag grn2.Dg];
numGenes = size (grn1.Tgg, 1);
numRegs  = size (grn1.Tgg, 2);
rowNames = opts.geneNames(1:numGenes);
colNames = [opts.geneNames(1:numRegs)  'h' 'R' '\lambda' 'D'];
numRows = numel (rowNames);
numCols = numel (colNames);



%set(gca,'FontSize',30);
set(0,'DefaultAxesFontSize',12);

hfig = figure(1);
set (hfig, 'pos',[0 800 800 600], 'Toolbar','None','MenuBar','None');

myColormap = .2*[5 0 0;5 1 1;5 2 2;5 3 3;5 4 4;5 5 5;4 5 4;3 5 3;2 5 2;1 5 1;0 5 0];
colormap (myColormap);

subplot (2,1,1);
imagesc (mat1, .05*[-1 1]);
colorbar;
%axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:numCols],'xticklabel', colNames);
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
title ('grn1: T_{gf}, h_g, R_g, \lambda_g, D_g');
for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat1(r, c) ), 'HorizontalAlignment','center','FontSize',14);
    end
end

subplot (2,1,2);
imagesc (mat2, .05*[-1 1]);
colorbar;
%axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:numCols],'xticklabel', colNames);
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
title ('grn2: T_{gf}, h_g, R_g, \lambda_g, D_g');
for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat2(r, c) ), 'HorizontalAlignment','center','FontSize',14);
    end
end
%print ('fig3', '-dpng');




return;


%======== T AND TInferred VERSUS 32 GENE INDICES
figure (4); clf; hold on;
title ( sprintf ('T-matrix elements') );
xlabel ('Element no.');
ylabel ('Element');
xlim ([0 numel(Tgg)])
ylim ([-100 200]);
grid on;
plot (Tgg(:), 'bx', 'LineWidth', 3);
plot (TggInferred(:), 'ro', 'LineWidth', 3);
legend ('Actual', 'Inferred');
x = linspace (-1000,1000); % Manually draw axes
y = linspace (0,0);
plot(x,y,'k-');
print ('fig4', '-dpng');

%======== COVARIANCE PLOT
figure (5); clf; hold on;
title ( sprintf ('Covariance of T vs Tinferred') );
xlabel ('T_actual');
ylabel ('T_inferred');
xlim ([-100 200]);
ylim ([-100 200]);
grid on;
plot (Tgg(:), TggInferred(:), 'ro', 'LineWidth', 3);
x = linspace (-1000,1000); % Manually draw axes
y = linspace (0,0);
plot(x,y,'k-');
plot(y,x,'k-');
print ('fig5', '-dpng');
end

