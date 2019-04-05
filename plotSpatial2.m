% Manu 11/24/18: Attempt to conform to Yen Lee's software design
% principles, but am not sure.

function [] = plotSpatial2 (opts, xntg, tt, g)
  numNuclei     = size(xntg, 1);
  numTimepoints = size(xntg, 2);
  
  nuclei = [1:1:numNuclei];

  %======= Resistor color codes for 1234567890 (9=white replaced by cyan)
  colorList = .2 * [3 2 1; 5 0 0; 5 3 0; 4 4 0; 0 5 0; 0 0 5; 5 0 5; 4 4 4; 0 5 5; 0 0 0];
  
  set(gcf,'DefaultAxesColorOrder',colorList)  
  plot (nuclei, xntg(:,:,g), 'LineWidth', 2);
  legend (sprintfc('t=%0.1f',tt), 'Position', [.8 .3 .1 .5]);
  xlabel ('Nucleus n');
  ylabel ('Expression x');
  title (sprintf ('Gene %s', opts.geneNames{g}));
end

