% this is a script where you can prepare code before you run it
addpath ('~/lib/glassmodels/subaxis');
close all;

%==== PLOT GRNs in one figure file using subaxis

hfig = figure();
set (gcf, 'Units', 'centimeters', 'Position',[0 0 16 14], 'Toolbar','None','MenuBar','None');

set(0,'DefaultAxesFontSize',10);
myhsv = [[ones(128,1); 0.4*ones(128,1)] 0.5*[[1:-1/127:0]'; [0:1/127:1]'] 0.9*ones(256,1)];
myColormap = hsv2rgb(myhsv);
%myColormap = .2*[5 2 2;5 2.5 2.5; 5 3 3;5 3.5 3.5;5 4 4;5 5 5;4 5 4;3 5 3;2 5 2;1 5 1;0 5 0];
colormap (myColormap);

geneNames = {'Hb', 'Kr', 'Gt', 'Kni', 'Bcd', 'Cad', 'Tll'};

%======== The CBI parameters after refinement
grntoplot = grnREF_sig;

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
rowNames = geneNames(1:numGenes);
colNames = [geneNames(1:numRegs)  'h' 'R' '\lambda' 'D'];
numRows = numel (rowNames);
numCols = numel (colNames);

subaxis (2,1,1,1, 'Spacing',0.02,'Padding',0,'ML',.05,'MR',.025,'MT',0.0,'MB',0.0);

matcolors = [grntoplot.Tgg grntoplot.hg zeros(numGenes,3)];
imagesc (matcolors, 0.2*[-1 1]);
axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:numCols],'xticklabel', colNames);
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
set(gca,'XAxisLocation', 'top');
for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat(r, c) ), 'HorizontalAlignment','center','FontSize',10);
    end
end



set (gca, 'FontSize', 10);
set (gca, 'FontWeight', 'bold');
set (gca, 'TickLength', [0 0.025]);
text(0.0,0.25,'A', 'FontSize', 12, 'FontWeight', 'bold');
%======== The SA parameters 
grntoplot = grnSA;

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
rowNames = geneNames(1:numGenes);
colNames = [geneNames(1:numRegs)  'h' 'R' '\lambda' 'D'];
numRows = numel (rowNames);
numCols = numel (colNames);

subaxis (2,1,1,2, 'Spacing',0.02,'Padding',0,'ML',.05,'MR',.025,'MT',0.0,'MB',0.0);

matcolors = [grntoplot.Tgg grntoplot.hg zeros(numGenes,3)];
imagesc (matcolors, .2*[-1 1]);
colorbar('Location','SouthOutside');
axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:numCols],'xticklabel', {});
set(gca,'ytick',[1:numRows],'yticklabel', rowNames);
for r=1:numRows
    for c=1:numCols
        text (c, r, sprintf ('%0.2f', mat(r, c) ), 'HorizontalAlignment','center','FontSize',10);
    end
end



set (gca, 'FontSize', 10);
set (gca, 'FontWeight', 'bold');
set (gca, 'TickLength', [0 0.025]);
text(0.0,0.25,'B', 'FontSize', 12, 'FontWeight', 'bold');


printpdf (gcf, 'grns-DrosCBISA.pdf');
return;

%=========== plot heatmaps of data and gap gene fits 
%==== PLOT xntgEXPT
genesymbols = {'Hunchback', 'Kr\"uppel', 'Giant', 'Knirps'};
hfig = figure ('Units', 'centimeters');
set (hfig, 'Position',[0 0 16 16], 'Toolbar','None','MenuBar','None');
for g=1:numGenes
    subaxis (5,numGenes,g,1, 'Spacing',0.02,'Padding',0,'ML',.1,'MR',.075,'MT',0.04,'MB',0.07);
    plotExpressionHeatmap (opts,xntgEXPT,tt,g) ;

    ah = gca;
    if g == 1
        set(ah, 'YTick', [2:2:8], 'YTickLabel', [3.1:13.5:43.6]); 
        ylabel('$t$ (min)', 'Interpreter', 'latex');
        text(-7,0,'A', 'FontSize', 12, 'FontWeight', 'bold');
    else
        set(ah, 'YTick', [2:2:8], 'YTickLabel', {});
    end
    set(ah, 'XTick', [10:10:50], 'XTickLabel', {});
    title(genesymbols{g}, 'Interpreter', 'latex');

    if g == 4
        oldpos = get(ah, 'Position');
        cbh = colorbar(ah, 'Location', 'manual', ...
            'Position', [oldpos(1)+oldpos(3)+0.01 oldpos(2) 0.01 oldpos(4)]);
    end    

end

%==== PLOT yntgEXPT (=diagnostics.yntg)
for g=1:numGenes
    subaxis (5,numGenes,g,2, 'Spacing',0.02,'Padding',0,'ML',.1,'MR',.075);
    plotOnoffstateHeatmap (opts,diagnostics.yntg,tt,g) ;

    ah = gca;
    if g == 1
        set(ah, 'YTick', [2:2:8], 'YTickLabel', [3.1:13.5:43.6]); 
        ylabel('$t$ (min)', 'Interpreter', 'latex');
        text(-7,0,'B', 'FontSize', 12, 'FontWeight', 'bold');
    else
        set(ah, 'YTick', [2:2:8], 'YTickLabel', {});
    end
    set(ah, 'XTick', [10:10:50], 'XTickLabel', {});

end

%==== PLOT yntgCBI (i.e., sgn (sum_g T_gf x_f + h_g) for xntgCBI)
numNuclei = size(xntgEXPT,1);
numTimepoints = size(xntgEXPT,2);
yntgMANU = zeros(numNuclei,numTimepoints,numGenes);
for n=1:numNuclei
    for t=1:numTimepoints
        yntgMANU(n,t,:) = sign (grnCBI.Tgg(:,:) * squeeze(xntgEXPT(n,t,:)) + grnCBI.hg(:) );
    end
end

for g=1:numGenes
    subaxis (5,numGenes,g,3, 'Spacing',0.02,'Padding',0,'ML',.1,'MR',.075);
    plotOnoffstateHeatmap (opts,yntgMANU,tt,g) ;

    ah = gca;
    if g == 1
        set(ah, 'YTick', [2:2:8], 'YTickLabel', [3.1:13.5:43.6]); 
        ylabel('$t$ (min)', 'Interpreter', 'latex');
        text(-7,0,'C', 'FontSize', 12, 'FontWeight', 'bold');
    else
        set(ah, 'YTick', [2:2:8], 'YTickLabel', {});
    end
    set(ah, 'XTick', [10:10:50], 'XTickLabel', {});



end

%==== PLOT xntgREF_sig - prediction of refined CBI model
for g=1:numGenes
    subaxis (5,numGenes,g,4, 'Spacing',0.02,'Padding',0,'ML',.1,'MR',.075);
    plotExpressionHeatmap (opts,xntgREF_sig,tt,g) ;

    ah = gca;
    if g == 1
        set(ah, 'YTick', [2:2:8], 'YTickLabel', [3.1:13.5:43.6]); 
        ylabel('$t$ (min)', 'Interpreter', 'latex');
        text(-7,0,'D', 'FontSize', 12, 'FontWeight', 'bold');
    else
        set(ah, 'YTick', [2:2:8], 'YTickLabel', {});
    end
    set(ah, 'XTick', [10:10:50], 'XTickLabel', {});

end

%==== PLOT xntgSA - prediction of simulated annealing model
for g=1:numGenes
    subaxis (5,numGenes,g,5, 'Spacing',0.02,'Padding',0,'ML',.1,'MR',.075);
    plotExpressionHeatmap (opts,xntgSA,tt,g) ;

    ah = gca;
    if g == 1
        set(ah, 'YTick', [2:2:8], 'YTickLabel', [3.1:13.5:43.6]); 
        ylabel('$t$ (min)', 'Interpreter', 'latex');
        text(-7,0,'E', 'FontSize', 12, 'FontWeight', 'bold');
    else
        set(ah, 'YTick', [2:2:8], 'YTickLabel', {});
    end
    set(ah, 'XTick', [10:10:50], 'XTickLabel', [45:10:85]);
    xlabel('$n$ (\% EL)', 'Interpreter', 'latex');

end
printpdf (gcf, 'heatmaps-DrosCBISA.pdf');
return;




%==== PLOT xntgCBI
hfig = figure ();
set (hfig, 'pos',[0 500 800 200], 'Toolbar','None','MenuBar','None');
for g=1:numGenes
    subaxis (1,numGenes,g, 'Spacing',.025,'Padding',0,'ML',.025,'MR',.025);
    plotExpressionHeatmap (opts,xntgCBI,tt,g) ;
end
printpdf (gcf, 'heatmap-xntgCBI.pdf');


%==== PLOT xntgEXPT
hfig = figure ();
set (hfig, 'pos',[0 500 800 200], 'Toolbar','None','MenuBar','None');
for g=1:numGenes
    subaxis (1,numGenes,g, 'Spacing',.025,'Padding',0,'ML',.025,'MR',.025);
    plotExpressionHeatmap (opts,xntgEXPT,tt,g) ;
end
printpdf (gcf, 'heatmap-xntgEXPT.pdf');

%==== PLOT yntgEXPT (=diagnostics.yntg)
hfig = figure ();
set (hfig, 'pos',[0 500 800 200], 'Toolbar','None','MenuBar','None');
for g=1:numGenes
    subaxis (1,numGenes,g, 'Spacing',.025,'Padding',0,'ML',.025,'MR',.025);
    plotOnoffstateHeatmap (opts,diagnostics.yntg,tt,g) ;
end
printpdf (gcf, 'heatmap-yntgEXPT.pdf');

%==== PLOT yntgCBI (i.e., sgn (sum_g T_gf x_f + h_g) for xntgCBI)
numNuclei = size(xntgEXPT,1);
numTimepoints = size(xntgEXPT,2);
yntgMANU = zeros(numNuclei,numTimepoints,numGenes);
for n=1:numNuclei
    for t=1:numTimepoints
        yntgMANU(n,t,:) = sign (grnCBI.Tgg(:,:) * squeeze(xntgEXPT(n,t,:)) + grnCBI.hg(:) );
    end
end
hfig = figure ();
set (hfig, 'pos',[0 500 800 200], 'Toolbar','None','MenuBar','None');
for g=1:numGenes
    subaxis (1,numGenes,g, 'Spacing',.025,'Padding',0,'ML',.025,'MR',.025);
    plotOnoffstateHeatmap (opts,yntgMANU,tt,g) ;
end
printpdf (gcf, 'heatmap-yntgMANU.pdf');
return;




%==== PLOT GRNs

%set(gca,'FontSize',14);
set(0,'DefaultAxesFontSize',12);

plotGRN (opts, grnSA)  ; title ('grnSA:  T_{gf}, h_g, R_g, \lambda_g, D_g');
set (gcf, 'pos',[0 430 600 200], 'Toolbar','None','MenuBar','None');
set (gca, 'FontSize', 12);
set (gca, 'FontWeight', 'bold');
set (gca, 'TickLength', [0 0.025]);
print ('grnSA', '-dpng', '-r300');

print ('grnSA', '-dpdf', '-r300');


plotGRN (opts, grnCBI) ; title ('grnCBI: T_{gf}, h_g, R_g, \lambda_g, D_g');
set (gcf, 'pos',[0 200 600 200], 'Toolbar','None','MenuBar','None');
set (gca, 'FontSize', 12);
set (gca, 'FontWeight', 'bold');
set (gca, 'TickLength', [0 0.025]);
print ('grnCBI-SLOPE', '-dpng', '-r300');

return;




%----- plot profile evolution
for g=1:numGenes
    hfig = figure(1);
    set (hfig, 'pos',[0 500 400 300], 'Toolbar','None','MenuBar','None');
    plotSpatial2 (opts, xntgEXPT, tt, g);
    set (gca, 'FontSize', 12);
    set (gca, 'FontWeight', 'bold');
    print (sprintf('spatial-xntgEXPT%d',g), '-dpng', '-r300');
    
    
    hfig = figure(2);
    set (hfig, 'pos',[500 500 400 300], 'Toolbar','None','MenuBar','None');
    plotSpatial2 (opts, xntgCBI, tt, g);
    set (gca, 'FontSize', 12);
    set (gca, 'FontWeight', 'bold');
    print (sprintf('spatial-xntgCBI-KINK%d',g), '-dpng', '-r300');
    
    pause
end
return;


%=========== I NEED THIS FUNCTION.
% THE BUILTIN MATLAB FUNCTIONS DON'T WORK PROPERLY.
%print ('heatmap-xntgEXPT', '-bestfit', '-dpdf', '-r300');
%saveas (gcf, 'heatmap-xntgEXPT.pdf');
%
function printpdf(h,outfilename)
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);
end
