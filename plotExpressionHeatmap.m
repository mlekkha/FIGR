function [] = plotExpressionHeatmap (opts,xntg,tt,g)

%======== HEATMAP OF GENE EXPRESSION TRAJECTORIES x_ntg
%'pos',[600 400 400 400], 
set (gcf, 'Toolbar','None','MenuBar','None');
imagesc(xntg(:,:,g)', [0 225]);
ah = gca;
colormap(ah, 'jet');

%title (sprintf ('Expression x(n,t) of gene %s', opts.geneNames{g}), 'FontSize', 16);
%xlabel ('Nucleus n', 'FontSize', 16); %AP position (% EL)
%ylabel ('Time class t', 'FontSize', 16);

%print (sprintf('heatmap-xntg%d',g), '-dpng', '-r300');
end

