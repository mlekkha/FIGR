function [] = plotExpression3DComparison (opts,xntg1,xntg2,tt,g)

numNuclei = size(xntg1,1);

%======== 3D PLOT OF GENE EXPRESSION TRAJECTORIES x_ntg
hfig = figure(1);
set (hfig, 'pos',[0 50 600 600]);
clf; hold on;

color1 = [0. 0. 1.];
color2 = [0. 1. 0.];
for n = 1:numNuclei
    %mySpline = csaps (tt, xntg1 (n,:,g), opts.splinesmoothing);
    nArray = n + 0*tt;
    plot3 (tt, nArray, xntg1 (n,:,g), 'Color', color1, 'LineWidth', 1);
    plot3 (tt, nArray, xntg2 (n,:,g), 'Color', color2, 'LineWidth', 2);
    scatter3 (tt, nArray, xntg1(n,:,g), 40,'o','MarkerEdgeColor', color1);
    scatter3 (tt, nArray, xntg2(n,:,g), 40,'o','MarkerEdgeColor', color2);
    
    %plot3 (tt, nArray, fnval(mySpline,tt), 'Color', [.6 .6 1.]);
    %scatter3 (tt, nArray,xntg(n,:,g), 40,'o','MarkerEdgeColor', [1 0 0]);
end
view (30, 45);
title (sprintf ('Expression x(n,t) of gene %s', opts.geneNames{g}), 'FontSize', 16);
xlabel ('Time t');
ylabel ('Nucleus n');
zlabel ('Expression x');
%print (sprintf('heatmap-xntg%d',g), '-dpng', '-r300');



figure(2); clf; hold on;
disp (size( squeeze(xntg1 (:,:,g))));
nn = (1:numNuclei)';
[nmesh, tmesh] = meshgrid (tt, nn);
disp (size(nmesh));
disp (size(tmesh));
%colormap (copper);
colormap (winter);
surf (nmesh, tmesh, squeeze(xntg1 (:,:,g)), 0+0*nmesh);
surf (nmesh, tmesh, squeeze(xntg2 (:,:,g)), 1+0*nmesh);
%shading interp
end

