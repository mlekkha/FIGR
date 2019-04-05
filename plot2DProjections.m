function [hfig] = plot2DProjections (opts,xntg,yntg,tt,g,g1,g2)

% Plot data points for gene g, classified as "on" or "off",
% in (g1,g2) space

numNuclei     = size (xntg,1);
numTimepoints = size (xntg,2);
numRegulators = size (xntg,3);
xkg = reshape (xntg, [numNuclei*numTimepoints numRegulators]);
ykg = reshape (yntg, [numNuclei*numTimepoints numRegulators]);
yk = ykg(:,g);

%======== Plot (a la nColScat2)
% move these to separate func and call from main()
negs = find (yk<=0);  % note that yk has already been clamped to 0
poss = find (yk>0);
hfig = figure();
set (hfig,'pos',[900 550 500 500],'Toolbar','None','MenuBar','None');
clf; hold on;

axis image;
xlim ([0 250]); ylim ([0 250]);
title (sprintf ('On/off data points for gene %s', opts.geneNames{g}), 'FontSize', 16);
xlabel (sprintf ('Gene %s', opts.geneNames{g1}), 'FontSize', 16);
ylabel (sprintf ('Gene %s', opts.geneNames{g2}), 'FontSize', 16);
scatter (xkg(negs,g1), xkg(negs,g2), 'o',  'MarkerEdgeColor', [0 .5 0]);
scatter (xkg(poss,g1), xkg(poss,g2), 'r*');
end

