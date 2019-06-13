function [] = plotDiscrepancies (runParams, debug)

%datadir = '~yenlee/genecircuits/';

%datadir = '';
datadir = '~/code/glassmodels/';

path (path(), '~/code/glassmodels/distributionPlot');

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontWeight','Normal');  % 'Bold'
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');

% Make figs larger by a factor of 2?
figure(1); clf;
set(gcf,'defaultAxesFontSize',12);
set (gcf, 'Units', 'centimeters', 'Position',[0 0 18 6],'Toolbar','None','MenuBar','None');

%========== ALL QUANTITIES, 2 genes, 10000 tmats, 10 nuclei ======================
% clf;
% discrep = dlmread ([datadir '2genes_10000tmats_10nuclei/discrep.dat']);
% distributionPlot (discrep, ...
%     'xnames', {'$\delta_T$', '$\delta_h$', '$\delta_R$', '$\delta_\lambda$'}, ...
%     'histOpt', 0, 'divFactor', 500, ...
%     'color', [.5 .5 .8], ...
%     'showMM', 6);
% set (gca, 'TickLabelInterpreter', 'latex');
% xlabel ('Quantity');
% ylabel ('Discrepancies');
% xlim ([0.4 4.6]); ylim ([0 0.4]);
% printpdf (gcf, [figdir '2genes_10000tmats_10nuclei_ThRl.pdf']);


%========== ALL QUANTITIES, 2 genes, 10000 tmats, 100 nuclei ======================
% DISCREPANCIES "discrep" ARE CALCULATED HERE, ON-THE-FLY
subplot(1,2,1);
grnTOY = dlmread ([datadir '2genes_10000tmats_100nuclei/grnTOY.dat']);
grnCBI = dlmread ([datadir '2genes_10000tmats_100nuclei/grnCBI.dat']);
%grnTOY = dlmread ([datadir '2genes_100tmats_30nuclei/grnTOY.dat']);
%grnCBI = dlmread ([datadir '2genes_100tmats_30nuclei/grnCBI.dat']);
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

%printpdf (gcf, [figdir '2genes_10000tmats_100nuclei_ThRl.pdf']);


%========== deltaT, 2 genes, 10000 tmats, 10/100/1000 nuclei ======================
% Only consider discrepancies in T; ignore h, R, lambda.
% It's best to accumulate the distributionPlots one at a time.
subplot(1,2,2);
tmp1 = dlmread ([datadir '2genes_10000tmats_10nuclei/grnTOY.dat']);  
tmp2 = dlmread ([datadir '2genes_10000tmats_10nuclei/grnCBI.dat']);  
dat1 = vecnorm(tmp1(:,1:2) - tmp2(:,1:2), 2, 2);

tmp1 = dlmread ([datadir '2genes_10000tmats_30nuclei/grnTOY.dat']);  
tmp2 = dlmread ([datadir '2genes_10000tmats_30nuclei/grnCBI.dat']);  
dat2 = vecnorm(tmp1(:,1:2) - tmp2(:,1:2), 2, 2);

tmp1 = dlmread ([datadir '2genes_10000tmats_100nuclei/grnTOY.dat']);  
tmp2 = dlmread ([datadir '2genes_10000tmats_100nuclei/grnCBI.dat']);  
dat3 = vecnorm(tmp1(:,1:2) - tmp2(:,1:2), 2, 2);

%tmp = dlmread ([datadir '4genes_10000tmats_10nuclei/discrep.dat']);  dat4 = tmp(:,1);
%tmp1 = dlmread ([datadir '4genes_10000tmats_10nuclei/grnTOY.dat']);  
%tmp2 = dlmread ([datadir '4genes_10000tmats_10nuclei/grnCBI.dat']);  
%dat4 = vecnorm(tmp1(:,1:4) - tmp2(:,1:4), 2, 2);

%tmp = dlmread ([datadir '4genes_1600tmats_100nuclei/discrep.dat']);  dat4 = tmp(:,1);
tmp1 = dlmread ([datadir '4genes_1600tmats_100nuclei/grnTOY.dat']);  
tmp2 = dlmread ([datadir '4genes_1600tmats_100nuclei/grnCBI.dat']);  
dat4 = vecnorm(tmp1(:,1:4) - tmp2(:,1:4), 2, 2);

%tmp1 = dlmread ([datadir '6genes_10000tmats_10nuclei/grnTOY.dat']);  
%tmp2 = dlmread ([datadir '6genes_10000tmats_10nuclei/grnCBI.dat']);  
%dat6 = vecnorm(tmp1(:,1:6) - tmp2(:,1:6), 2, 2);

%tmp = dlmread ([datadir '6genes_1600tmats_100nuclei/discrep.dat']);  dat5 = tmp(:,1);
tmp1 = dlmread ([datadir '6genes_1600tmats_100nuclei/grnTOY.dat']);  
tmp2 = dlmread ([datadir '6genes_1600tmats_100nuclei/grnCBI.dat']);  
dat5 = vecnorm(tmp1(:,1:6) - tmp2(:,1:6), 2, 2);

distributionPlot (dat1, 'xValues',1, 'histOpt',0, 'divFactor',500, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat2, 'xValues',2, 'histOpt',0, 'divFactor',500, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat3, 'xValues',3, 'histOpt',0, 'divFactor',500, 'color', [.5 .5 .8], 'showMM',6);
distributionPlot (dat4, 'xValues',4, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);
distributionPlot (dat5, 'xValues',5, 'histOpt',0, 'divFactor',100, 'color',[.5 .5 .8], 'showMM',6);

ylabel ('\delta_T');
xlim ([0.5 5.+.5]); ylim ([-0.05 1.0]);
set (gca, 'TickLabelInterpreter', 'latex');
xticks ([1 2 3 4 5]);
xticklabels ({ ...
    '$\matrix{N=10 \cr G=2}$'  , ...
    '$\matrix{N=30 \cr G=2}$'  , ...
    '$\matrix{N=100 \cr G=2}$'  , ...
    '$\matrix{N=100 \cr G=4}$'  , ...
    '$\matrix{N=100 \cr G=6}$'   });

text(-0.16,1,'B', 'FontSize', 12, 'FontWeight', 'bold');
ah2 = gca;

set(ah1, 'Units', 'centimeters');
set(ah2, 'Units', 'centimeters');
set(ah1, 'Position', [1.5 1.5 6.67 4]);
set(ah2, 'Position', [9.67 1.5 8.33 4]);

export_fig -transparent -nocrop ~manu/docs/manuscripts/GRN_ms/SyntheticDataTestsByRow.pdf

%text (1, .2, '$\matrix{ N=10 \cr sdf}$', 'FontSize', 14, 'interpreter', 'latex');

%set(get(gca,'YLabel'),'Rotation',0);
%ylh = get(gca, 'ylabel'); ylp = get(ylh, 'Position');
%set (ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
%printpdf (gcf, [figdir 'XXXgenes_10000tmats_XXXnuclei_T.pdf']);











return;

discrep2 = dlmread ([datadir '2genes_10000tmats_10nuclei/discrep.dat']);
discrep4 = dlmread ([datadir '~/genecircuits/4genes_1000tmats_100nuclei/discrep.dat']);
discrep6 = dlmread ([datadir '~/genecircuits/6genes_1000tmats_100nuclei/discrep.dat']);
discrep = [discrep2 ; discrep4 ; discrep6];
discrep (isnan(discrep)) = -0.2;   % Represent failed inferences by this value of delta
groups  = [repmat(2,size(discrep2)) ; repmat(4,size(discrep4)) ;  repmat(6,size(discrep6)) ];

figure(3); clf; hold on;
distributionPlot (discrep,'groups',groups,'colormap',[.2 .8 .8],'showMM',6);
text (1, 1.3, '1000 T-matrices, 100 trajectories each, 41 timepoints', 'FontSize', 14);
xlabel('Number of genes');
ylabel('Discrepancy delta');
xlim ([1.0 7.0]);
ylim ([-0.3 1.5]);
%title ('Quality of Parameter Inference Decreases with Number of Genes');
%text (0.2, -0.1, 'Failures ->      ', 'FontSize', 14);
myYTick = get (gca, 'YTick');
myYTickLabel = get (gca, 'YTickLabel');
myYTick = [-0.2  myYTick];
myYTickLabel = char ('Fail' , myYTickLabel);
set (gca, 'YTick', myYTick);
set (gca, 'YTickLabel', myYTickLabel);
print ('discrep-vs-numgenes','-dpng', '-r300'); % 300 DPI

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

