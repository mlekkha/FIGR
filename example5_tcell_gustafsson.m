    %% FIGR example for hematopoietic FDCIP-mix data from May et al's experiment.

clc;

%======== READ NCBI FILES? =======
%[xntgEXPT tt nucleusNames geneNames]    ...
%= loadNCBIFiles("eryneu_xntg.txt", "eryneu_tt.txt");

% CONVERT TO MDA FORMAT?
%writematrix (geneNames, "eryneu_gnames.txt");
%saveMDA ("eryneu_xntg.mda", xntgEXPT);  % this is really just optional
%saveMDA ("eryneu_tt.mda", tt);          % this is really just optional

% CONVERT TO PLAIN TEXT?  (NOT GOOD BECAUSE 3D ARRAYS NOT SUPPORTED)
%writematrix (tt, 'eryneu_tt.mdatxt', 'FileType', 'text');
%writematrix (xntgEXPT, 'eryneu_xntg.mdatxt', 'FileType', 'text');

%======== READ EXPERIMENTAL TRAJECTORIES xntg(:,1,:) AND TIMEPOINTS tt =======
[xntgEXPT,tt,nucleusNames,geneNames] = ...
                loadNCBIFiles ('tcell-gustafsson-xntg.txt', ...
                               'tcell-gustafsson-tt.txt');

%======== READ VALUES OF OPTIONS p, v, and x THAT HAVE BEEN TUNED BY USER  =======
% INDEX ORDER IS n, g, o (nucleus, gene, option index)
pvxOpts_ngo = loadMDA ('tcell_options_06_02_21_1.mda');

%======== Define global structs for options and ODE options
global opts;
global ODEopts;
global optimopts;

%======== Start timer
tic;

%======== SET OPTIONS (see README.md for description) ========
% NOTE: slopethresh, exprthresh, splinesmoothing are no longer used!
% NOTE: Rld_tsafety also should be removed.
% NOTE: These pars are now supplied via  pvxOpts_ngo
opts = struct(  'debug', 0, ...
    'slopethresh', NaN, ...
    'exprthresh', NaN, ...
    'splinesmoothing', NaN, ...
    'Rld_tsafety', 0, ...       % should eventually ged rid
    'spatialsmoothing', 0.5, ...
    'minborder_expr_ratio', 0.01, ...
    'Rld_method', 'slope_nodiff', ...
    'synthesisfunction', 'synthesis_sigmoid_sqrt', ...
    'ODEAbsTol', 1e-3, ...
    'ODEsolver', 'ode45', ...
    'pvxOpts_ngo', pvxOpts_ngo, ...
    'lambda', 0.5, ...
    'lm', 'FIGRlogReg'); % glmfit|FIGRlogReg|lassoglm

%======== set the ODE options
ODEopts = odeset('AbsTol', opts.ODEAbsTol);

numGenes = 12;  % have to hardwire
numNuclei = size (xntgEXPT,1);

%======== INFER GRN PARAMETERS grnFIGR
[grnFIGR, diagnostics] = infer (opts, xntgEXPT, tt, numGenes);
yntgEXPT = diagnostics.yntg;
saveMDA ("tcell-gustafsson-yntg.mda", yntgEXPT);


%======== REFINE GRN BY NELDER-MEAD ========
xntgFLAT = reshape(xntgEXPT, 12, 13);

% set optimization options for refinement
packed_paramvec = packParams(grnFIGR, numGenes);
optimopts = optimset( 'Display', 'Iter', ...
                      'MaxFunEvals', 200*length(packed_paramvec), ...
                      'MaxIter', 200*length(packed_paramvec));
                        
% Refine parameters
grnREF = refineFIGRParams(grnFIGR, xntgFLAT, tt);


%======== Print status
disp ('Nelder-Mead refinement complete... ');

%======== Print time taken
fprintf(1, 'Total time elapsed: %f\n', toc);

%======== RECOMPUTE TRAJECTORIES
[xntgRECAL] = computeTrajs (opts, grnREF, xntgEXPT, tt);
  

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

% save T matrix
filename = 'tcell-gustafsson-Tmatrix.txt';
saveTMatrix (filename, grnREF.Tgg, geneNames);

% save the plots
exportgraphics(gcf,'tcell-gustafsson-fits.pdf','BackgroundColor','none');


function saveTMatrix  (filename, a, geneNames)
delim = {'\t', '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n', 'n', 'n', 'n'};
fmtstr = '%.2f';
a = a';
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
