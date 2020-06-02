%% FIGR example for hematopoietic FDCIP-mix data from May et al's experiment.

clc;
 
%======== READ EXPERIMENTAL TRAJECTORIES xntg(:,1,:) AND TIMEPOINTS tt =======
[xntgEXPT tt geneNames] = readGeneExprFiles();

%======== READ PARAMETERS p, v, and x DEFINED BY USER =======
pvxParams = readArray('params.txt');

%======== Define global structs for options and ODE options 
global opts;
global ODEopts;
global optimopts;

%======== SET OPTIONS (see README.md for description) ========
opts = struct(  'debug', 0, ...                 
                'slopethresh', 1.0, ...         
                'exprthresh', 0.5, ...
                'splinesmoothing', 0.01, ...
                'spatialsmoothing', 0.5, ...
                'minborder_expr_ratio', 0.01, ...
                'Rld_method', 'slope', ...
                'Rld_tsafety', 3, ...
                'synthesisfunction', 'synthesis_sigmoid_sqrt', ...
                'ODEAbsTol', 1e-3, ...
                'ODEsolver', 'ode45', ...
                'pvxParams', pvxParams, ...
                'lambda', 0.5, ...
                'lm', 'FIGRlogReg'); % glmfit, FIGRlogReg, lassoglm
            
%======== set the ODE options
ODEopts = odeset('AbsTol', opts.ODEAbsTol);

numGenes = 12; 

numNuclei = size (xntgEXPT,1);

%======== Start timer
tic;

%======== INFER GRN PARAMETERS grnFIGR
[grnFIGR, diagnostics] = infer (opts, xntgEXPT, tt, numGenes);
yntgEXPT = diagnostics.yntg;

  marker_size = 4;
   

  [xntgREF] = computeTrajs (opts, grnFIGR, xntgEXPT, tt);
  
  %[xntgREF] = computeTrajs (opts, grnREF, xntgEXPT, tt);
  close all;
  figure('Units', 'inches', 'Position', [0 0 8.5 5.75]);
  for g=1:numGenes
      for n=1:numNuclei
         h(g) = subplot(4,3,g);
              
         %set(a,'box','off','color','none')
         %b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
         if(n == 1)
            ph(1) = plot(tt, xntgREF(n,:,g), '-', 'MarkerFaceColor', 'r', 'color', 'r','LineWidth', 1.5);
         else
            ph(1) = plot(tt, xntgREF(n,:,g), '-', 'MarkerFaceColor', 'b', 'color', 'b','LineWidth', 1.5);
         end
         hold on;
         if(n == 1)
            ph(2) = plot(tt, xntgEXPT(n,:,g), 'o', 'MarkerFaceColor', 'r', 'color', 'r','markers', marker_size);
         else
            ph(2) = plot(tt, xntgEXPT(n,:,g), 'o', 'MarkerFaceColor', 'b', 'color', 'b','markers', marker_size);
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

filename = 'FIGR_Tmatrix.txt';
writeArray(filename, grnFIGR.Tgg, geneNames);

% writeArray() and readArray() provide a way to write arbitrary-dimensional
% arrays.
function writeArray (filename, a, geneNames)
delim = {'\t', '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n', 'n', 'n', 'n'};
fmtstr = '%.2f';
dmax = ndims(a);
cmaxd = size(a);
fid = fopen (filename, 'w');
for i=1:numel(geneNames)
    fprintf (fid, '%s', string(geneNames(i)));
    fprintf (fid, '\t');
end
fprintf (fid, '\n');
idx = ones (cmaxd);                 % start all indices at 1
for i=1:numel(a)
    fprintf (fid, fmtstr, a(i));      % print array element
    for d=1:dmax
        fprintf (fid, delim{d});    % print delimiter
        idx(d) = idx(d)+1;          % increment index
        if (idx(d) <= cmaxd(d)); break; end;
        idx(d) = 1;                 % restart index at 1
    end                             % consider index at next level
end
fclose (fid);
end


function [a] = readArray (filename)
fid = fopen (filename, 'r');
dmax = fscanf (fid, '%d');   % read number of dimensions
fgets (fid);                 % skip to next line
for d=1:dmax
    cmaxd(d) = fscanf (fid, '%d');  % read number of elements
    fgets (fid);                    % skip to next line
end
a = fscanf (fid, '%f');      % read in all elements
a = reshape (a, cmaxd);
fclose (fid);
end

function phasePlot (xntgEXPT, xntgREF, geneNames)

    n1 = 1;
    g1 = 1;
    n2 = 2;
    g2 = 2;
    plot(xntgEXPT(n1,:,g1), xntgEXPT(n2,:,g2));
    xGene = geneNames(g1);
    yGene = geneNames(g2);
    
    xlabel("df")
    ylabel("df")
    %plot3(xntgEXPT(1,:,1), xntgEXPT(2,:,2), xntgEXPT(2,:,3));

end