function [] = plotOnoffstateHeatmap (opts,yntg,tt,g)

%======== HEATMAP OF ON/OFF STATE y_ntg
imagesc(yntg(:,:,g)');
ah = gca;
colormap(ah, [[0.75:-0.01:0]', [0:0.01:0.75]', zeros(76,1)]);

%title (sprintf ('State y(n,t) of gene %s', opts.geneNames{g}), 'FontSize', 16);
%xlabel ('Nucleus n', 'FontSize', 16); %AP position (% EL)
%ylabel ('Time t (min)', 'FontSize', 16);
%print (sprintf('heatmap-yntg%d',g), '-dpng', '-r300');
end
