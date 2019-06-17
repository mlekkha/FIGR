path (path(), '~/code/glassmodels/distributionPlot');

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontWeight','Normal');  % 'Bold'
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');

close all;
figure(1); clf;
set(gcf,'defaultAxesFontSize',12);
%set (gcf, 'Units', 'centimeters', 'Position',[0 0 18 6],'Toolbar','None','MenuBar','None');
set (gcf, 'Units', 'centimeters', 'Position',[0 0 50 20],'Toolbar','None','MenuBar','None');

%========== ALL QUANTITIES, 2 genes, 10000 tmats, 100 nuclei ======================
% DISCREPANCIES "discrep" ARE CALCULATED HERE, ON-THE-FLY
subplot(1,2,1);
datadir = sprintf ('2genes_10000tmats_30nuclei/');
grnTOY = dlmread ([datadir '/grnTOY.dat']);
grnCBI = dlmread ([datadir '/grnCBI.dat']);
discrepT = vecnorm(grnTOY(:,1:2) - grnCBI(:,1:2), 2, 2);
discrep = [discrepT abs(grnTOY(:,3:5)-grnCBI(:,3:5))];
distributionPlot (discrep, ...
    'xnames', {'$\delta_T$', '$\delta_h$', '$\delta_R$', '$\delta_\lambda$'}, ...
    'histOpt', 0, 'divFactor', 500, ...
    'color', [.5 .5 .8], ...
    'showMM', 6);
set (gca, 'TickLabelInterpreter', 'latex');
xlabel ('$N=100, G=2$', 'Interpreter','latex', 'FontSize', 12);
ylabel ('Discrepancies', 'FontSize', 12);
xlim ([0.4 4.6]); ylim ([-0.015 0.3]);
text(-0.3,0.3,'A', 'FontSize', 12, 'FontWeight', 'bold');
ah1 = gca;



%========== deltaT, 2 genes, 10000 tmats, 10/100/1000 nuclei ======================
% Only consider discrepancies in T; ignore h, R, lambda.
% It's best to accumulate the distributionPlots one at a time.
subplot(1,2,2);

%dat1 = calcDiscrepancies (2, 10000, 10);
%dat2 = calcDiscrepancies (2, 10000, 30);
%dat6 = calcDiscrepancies (10, 100, 30);
dat1 = calcDiscrepancies (2, 10000, 100);
dat2 = calcDiscrepancies (4, 1600, 100);
dat3 = calcDiscrepancies (6, 1600, 100);
dat4 = calcDiscrepancies (15, 100, 100);
dat5 = calcDiscrepancies (20, 100, 100);
dat7 = calcDiscrepancies (40, 100, 100);
dat8 = calcDiscrepancies (50, 100, 100);

%return;

distributionPlot (dat1, 'xValues',1, 'histOpt',0, 'divFactor',500, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat2, 'xValues',2, 'histOpt',0, 'divFactor',500, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat3, 'xValues',3, 'histOpt',0, 'divFactor',500, 'color', [.5 .5 .8], 'showMM',6);
distributionPlot (dat4, 'xValues',4, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat5, 'xValues',5, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);
%distributionPlot (dat6, 'xValues',6, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat7, 'xValues',7, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat8, 'xValues',8, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);

xticklabels ({ ...
    %'$\matrix{N=10 \cr G=2}$'  , ...
    %'$\matrix{N=30 \cr G=2}$'  , ...
    '$\matrix{N=100 \cr G=2}$'  , ...
    '$\matrix{N=100 \cr G=4}$'  , ...
    '$\matrix{N=100 \cr G=6}$'  , ...
    '$\matrix{N=100 \cr G=15}$'  , ...
    '$\matrix{N=100 \cr G=20}$'  , ...
    '$\matrix{N=100 \cr G=??}$'  , ...
    '$\matrix{N=100 \cr G=40}$'  , ...
    '$\matrix{N=100 \cr G=50}$'  , ...
    ''   });
xticks (1:8);
set (gca, 'TickLabelInterpreter', 'latex');
ylabel ('\delta_T'); xlim ([0.5 9.+.5]); ylim ([-0.05 2.0]);
text(-0.16,1,'B', 'FontSize', 12, 'FontWeight', 'bold');
ah2 = gca;



%set(ah1, 'Units', 'centimeters');
%set(ah2, 'Units', 'centimeters');
%set(ah1, 'Position', [1.5 1.5 6.67 4]);
%set(ah2, 'Position', [9.67 1.5 8.33 4]);
%export_fig -transparent -nocrop ~manu/docs/manuscripts/GRN_ms/SyntheticDataTestsByRow.pdf


%======================= CALCULATE DISCREPANCIES FROM GRN PAR FILES ======================
function discreps = calcDiscrepancies (genes, tmats, nuclei)
datadir = sprintf ('%dgenes_%dtmats_%dnuclei/', genes, tmats, nuclei);
grnTOY = dlmread ([datadir '/grnTOY.dat']);
grnCBI = dlmread ([datadir '/grnCBI.dat']);
discreps = vecnorm(grnTOY(:,1:genes) - grnCBI(:,1:genes), 2, 2);
rows = numel(discreps);
nans = sum(isnan(discreps));
frac = nans/rows;
fprintf ('%d tmats \t%d nuclei\t%d genes:\t\t %d rows,\t%d NaNs,\t%f fracNaNs\n', tmats, nuclei, genes, rows, nans, frac);
end

%======================= WRITE TO PDF FILE ======================
function printpdf (h,outfilename)
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);
end

