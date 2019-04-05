function [] = plotExpression3D(opts,xntg,tt,yntg,g)

ynt = yntg (:,:,g);
numNuclei = size(xntg,1);

%======== 3D PLOT OF GENE EXPRESSION TRAJECTORIES x_ntg
hfig = figure(1);
set (hfig, 'pos',[0 50 600 600]);
clf; hold on;
for n = 1:numNuclei
    mySpline = csaps (tt, xntg(n,:,g), opts.splinesmoothing);
    nArray = n + 0*tt;
    negs = find (ynt(n,:) <= 0);
    poss = find (ynt(n,:) > 0);
    %plot3 (tt, nArray, xntg(n,:,g), 'Color', [.6 .6 1.]);
    plot3 (tt, nArray, fnval(mySpline,tt), 'Color', [.6 .6 1.]);
    scatter3 (tt(negs), nArray(negs),xntg(n,negs,g), 40,'o','MarkerEdgeColor', [1 0 0]);
    scatter3 (tt(poss),  nArray(poss),xntg(n,poss,g),40,'*','MarkerEdgeColor', [0 .5 0]);
end
view (30, 45);
title (sprintf ('Expression x(n,t) of gene %s', opts.geneNames{g}), 'FontSize', 16);
xlabel ('Time t');
ylabel ('Nucleus n');
zlabel ('Expression x');
%print (sprintf('heatmap-xntg%d',g), '-dpng', '-r300');
end

