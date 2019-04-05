% NOT UPDATED
% Manu 11/24/18: Attempt to conform to Yen Lee's software design
% principles, but am not sure.

function [] = plotSpatial (opts, xntgEXPT, xntgCBI, tt, numGenes)
  numNuclei     = size(xntgEXPT, 1);
  numTimepoints = numel (tt);
  numExternals  = size(xntgEXPT, 3) - numGenes;     

  nuclei = [1:1:numNuclei];

  %======= Resistor color codes for 1234567890
  colorList = .2 * [3 2 1; 5 0 0; 5 3 0; 4 4 0; 0 5 0; 0 0 5; 5 0 5; 4 4 4; 0 5 5; 0 0 0];

  for t=1:numTimepoints

    %-------- PLOT SIMULATED TRAJECTORIES OF HKGN 
    figure (1); clf; hold on;
    title ( sprintf ('Simulated concentrations of H/K/G/N, x(n,t=%g,i)', tt(t) ) );
    xlabel ('Nucleus n');
    ylabel ('Concentrations x');
    ylim ([0 250]);
    for i=1:numGenes
        theColor = colorList(i,:);
        plot (nuclei, squeeze(xntgCBI(:,t,i)), '-', ...
            'LineWidth', 3, 'Color', theColor, 'MarkerFaceColor', theColor);
    end
    legend (opts.geneNames{1:numGenes});
    %-------- ON SAME PLOT, COMPARE RAW TRAJECTOREIS OF HKGN 
    for i=1:numGenes
        theColor = colorList(i,:);
        plot (nuclei, squeeze(xntgEXPT(:,t,i)), '--x', ...
            'LineWidth', 1, 'Color', theColor, 'MarkerFaceColor', theColor);
    end
 
    %-------- PLOT GIVEN BCT
    figure (2); clf; hold on;
    title ( sprintf ('Measured concentrations of B/C/T, X(n,t=%g,i)', tt(t) ) );
    xlabel ('Nucleus n');
    ylabel ('Concentrations x');
    ylim ([0 250]);
    for i=numGenes+1:numGenes+numExternals
        theColor = colorList(i,:);
        plot (nuclei, squeeze(xntgEXPT(:,t,i)), '-', ...
            'LineWidth', 3, 'Color', theColor, 'MarkerFaceColor', theColor);
    end
    legend (opts.geneNames{numGenes+1:numGenes+numExternals});

    pause
  end

  fprintf ('=======================================\n\n\n\n');
end

