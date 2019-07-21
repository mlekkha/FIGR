% this is a script where you can prepare code before you run it
close all;


%=========== plot gap gene data, FIGR fits, and SA fits in 1D, time points
%are different panels
genesymbols = {'Hunchback', 'Kr\"uppel', 'Giant', 'Knirps'};
hfig = figure ('Units', 'centimeters');
set (hfig, 'Position',[0 0 16 16], 'Toolbar','None','MenuBar','None');
set(0,'DefaultAxesFontSize',8);
lcols = lines;
for tpt=1:8
    for g=1:numGenes
       
        subplot(8,numGenes,(tpt-1)*numGenes+g);

        hold on;
        pl1 = plot([35:92],xntgEXPT(:,tpt+1,g),'Color',[0.569 0.137 0.9], 'LineWidth', 1);  
        pl2 = plot([35:92],xntgSA(:,tpt+1,g),'Color',[0.9 0.569 0.137], 'LineWidth', 1);  
        pl3 = plot([35:92],xntgREF(:,tpt+1,g),'Color',[0.137 0.9 0.569], 'LineWidth', 1);  
        hold off;
        ylim([0,225]); 
        xlim([34,93]);

        ah = gca;
        if g == 1
            set(ah, 'YTick', [0:50:200], 'YTickLabel', {'0','','100','','200'}, 'FontSize', 8); 
            ylabel('$x$', 'Interpreter', 'latex');
            text(70,190,sprintf('$t = %.1f$', tt(tpt+1)), ...
                    'Interpreter', 'latex', 'FontSize', 8);
        else
            set(ah, 'YTick', [0:50:200], 'YTickLabel', {}, 'FontSize', 8);
        end

        if tpt == 8
            set(ah, 'XTick', [45:10:85], 'FontSize', 8); 
            xlabel('$n$ (\% EL)', 'Interpreter', 'latex');
        else
            set(ah, 'XTick', [45:10:85], 'XTickLabel', {}, 'FontSize', 8);
        end

        if tpt == 1
            title(genesymbols{g}, 'Interpreter', 'latex');
        end    

        if (tpt == 1) & (g == 2)

            legend(pl1, 'Data');
            legend('boxoff');

        end

        if (tpt == 1) & (g == 3)

            legend(pl2, 'SA');
            legend('boxoff');

        end

        if (tpt == 1) & (g == 4)

            legend(pl3, 'FIGR');
            legend('boxoff');

        end

        if (g == 3) & (tpt == 2) 

            text(82,2,'*');

        end

        if (g == 1) & (tpt == 7) 

            text(69,30,'$\triangleright$', 'Interpreter', 'latex');
            text(51,40,'$\triangleleft$', 'Interpreter', 'latex');

        end


    end
end

printpdf (gcf, 'lineplots-DrosCBISA.pdf');
return;




function printpdf(h,outfilename)
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);
end
