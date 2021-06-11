function vk = neyda(opts, xntgEXPT, grn, numGenes, numNuclei, tt, geneNames, yntgEXPT, xntgRECAL)

%numGenes = numel(grn.Rg);                  % G
numRegulators = size(grn.Tgg, 2);          % F
%numExternals = numRegulators - numGenes;   % F-G
numExternals = 0;

Tgf = grn.Tgg;   % INCLUDE ALL GENES
hg = grn.hg;
Rg = grn.Rg;
lambdag = grn.lambdag;
Dg = grn.Dg;


        numTimepoints = size(tt, 1);
        
        % accumulates the results from g function
        y = nan (numNuclei, numTimepoints, numGenes);
        
        % Synthesis
        for n=1:numNuclei %Number of Conditions

            for g=1:numGenes 
                
                for t=1:numTimepoints

                    if strcmp(opts.synthesisfunction, 'synthesis_sigmoid_sqrt')
                        xntg = reshape(xntgEXPT(n,t,:), [1, size(xntgEXPT,3)]); %reshape to concatonate into a vector
                        y(n,t,g) = synthesis_sigmoid_sqrt(Tgf(g,:) *  xntg' + hg(g)); 
                    elseif strcmp(opts.synthesisfunction, 'synthesis_heaviside')                   
                        xntg = reshape(xntgEXPT(n,t,:), [1, size(xntgEXPT,3)]); %reshape to concatonate into a vector
                        y(n,t,g) = synthesis_heaviside(Tgf(g,:) *  xntg' + hg(g)); 
                    else
                        error('Unknown synthesis function %s',...
                                    opts.synthesisfunction);
                    end
                
                end
   
            end
            
        end

%yntgEXPT(yntgEXPT == -1) = 0;
yntgEXPT(yntgEXPT == -1) = -2;
yntgEXPT(yntgEXPT == 1) = 2;
yntgEXPT(yntgEXPT == 0) = 3;
y(y == 0) = -1.5;
y(y == 1) = 1.5;

close all;
figure('Units', 'inches', 'Position', [0 0 8.5 5.75]);
marker_size = 4;
for g=1:numGenes
    for n=1:numNuclei
        h(g) = subplot(4,3,g);
        
        %set(a,'box','off','color','none')
        %b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
        if(n == 1)
           ph(1) = plot(tt, xntgRECAL(n,:,g), '-', 'MarkerFaceColor', 'r', 'color', 'r','LineWidth', 1.5);
        else
           ph(1) = plot(tt, xntgRECAL(n,:,g), '-', 'MarkerFaceColor', 'b', 'color', 'b','LineWidth', 1.5);
        end
        hold on;
        if(n == 1)
            ph(2) = plot(tt, xntgEXPT(n,:,g), 'o', 'MarkerFaceColor', 'r', 'color', 'r','markers', marker_size);
        
            ph(2) = plot(tt, y(n,:,g), '*', 'MarkerFaceColor', 'r', 'color', 'r','markers', 9);
        
            ph(2) = plot(tt, yntgEXPT(n,:,g), 's', 'MarkerFaceColor', 'r', 'color', 'r','markers', marker_size);
            
        else
            ph(2) = plot(tt, xntgEXPT(n,:,g), 'o', 'MarkerFaceColor', 'b', 'color', 'b','markers', marker_size);
            
            %ph(2) = plot(tt, y(n,:,g), 'o', 'MarkerFaceColor', 'b', 'color', 'b','markers', marker_size);
            
            ph(2) = plot(tt, y(n,:,g), '*', 'MarkerFaceColor', 'b', 'color', 'b','markers', 9);
        
            ph(2) = plot(tt, yntgEXPT(n,:,g), 's', 'MarkerFaceColor', 'b', 'color', 'b','markers', marker_size);
        end
        title(geneNames(g), 'FontSize', 10, 'fontweight', 'bold', 'interpreter','latex');
        
        if(g == 3)
            %legend([ph(2) ; ph(1)], {'Experimental data', 'Model'}, 'Location','northwest','Orientation','vertical', 'FontSize',8)
            %legend boxoff;
        end
        
        if(g == 10 || g == 11 || g == 12)
            % xlabel('$t$','FontSize', 12, 'interpreter', 'latex');
            xlabel('$t$','interpreter', 'latex');
        end
        if(g == 1 || g == 2)
            %set(gca,'xticklabel',{[]});
        end
        
        % axis([0 46.88 0 255]);
    end
end

        
        

        
        


end % from RateofChange()

function g = synthesis_sigmoid_sqrt(u)
g = 0.5*(u/sqrt(1.0+u*u) + 1.0);
end % from synthesis_sigmoid_sqrt()

function g = synthesis_heaviside(u)
if u >= 0.0 % Heaviside function,
    g = 1.0; %differs from inbuilt MATLAB one in that g(u) = 1
else         % instead of 0.5 if u = 0.
    g = 0.0;
end
end % from synthesis_heaviside()