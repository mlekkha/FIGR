clc;
fprintf (['\n'...
    '============================== Example 3 ==============================\n'...
    'Experimentally determined gene expression trajectories for seven genes \n'...
    '  (Hb,Kr,Gt,Kni,Bcd,Cad,Tll or HKGNBCT for short), together with timepoints,\n'...
    '  are read in from files fly_xntg.txt and fly_tt.txt.\n'...
    'FIGR is used to infer parameters of the gene regulatory network,\n'...
    '  assuming that HKGN regulate each other but BCT are upstream regulators.\n'...
    'This script visualizes the gene expression trajectories x(n,t,g) as heatmaps in (n,t) space\n'...
    '  where n=nucleus index and t=timepoint index.\n'...
    'The top panels show the experimentally measured data for H,K,G,N, \n'...
    '  whereas the bottom panels show the data re-simulated from the FIGR-inferred GRN.\n'...
    '\n'...
    'Although the re-simulated data do not match the experimental data perfectly,\n'...
    '  FIGR gives good starting parameters for refinement (e.g., Nelder-Mead), as in the paper.\n'...
    '=======================================================================\n'...
    '\n']);

% MANU TO DO:
% - CORRECT slopethresh AND exprthresh BELOW 
% - COULD ADD REFINEMENT EVENTUALLY
% - IS THERE A WAY TO SHARE ONE COLOR BAR FOR ALL HEATMAPS?


%======== SET OPTIONS FOR COMPUTETRAJS ========
opts.synthesisfunction = 'synthesis_sigmoid_sqrt';
opts.ODEAbsTol = 1e-2;                          % ODE solver tolerance
opts.ODEsolver = 'ode45';                       % 4th order Runge-Kutta
opts.debug = 0;                                 % verbosity level (0-3)
%======== SET OPTIONS FOR FIGR ========
opts.slopethresh = 1.0;  % YL 2018-11-6  -----need to see what Manu used last
opts.exprthresh = 100.0; % YL 2018-11-6
opts.splinesmoothing = 0.01;
opts.spatialsmoothing = 0.5;
opts.minborder_expr_ratio = 0.01;
opts.Rld_method = 'slope';
opts.Rld_tsafety = 3;
numGenes = 4;  % 4 genes are "internal"; 3 genes are upstream regulators

%======== READ EXPERIMENTAL TRAJECTORIES xntg(:,1,:) AND TIMEPOINTS tt =======
xntgEXPT = readArray ('fly_xntg.txt');  % xntgEXPT is now a 58x9x7 array
tt       = readArray ('fly_tt.txt');

%======== INFER GRN PARAMETERS grnFIGR
[grnFIGR, diagnostics] = infer (opts, xntgEXPT, tt, numGenes);
yntgEXPT = diagnostics.yntg;

%======== MANU MIGHT WANT TO PUT 'REFINE GRN BY NELDER-MEAD' HERE ========



%======== RECOMPUTE TRAJECTORIES xntg ========
[xntgFIGR] = computeTrajs (opts, grnFIGR, xntgEXPT, tt);

%======== DISPLAY grnSA AND grnFIGR FOR COMPARISON
disp ('GRN parameters inferred using FIGR: ');
disp ('grnFIGR = ');
disp (struct2table (grnFIGR));
disp ('grnFIGR.Tgg = ');
tmp = grnFIGR.Tgg;

opts.geneNames = {'H','K','G','N','B','C','T'};  % HARDWIRED
disp (array2table (tmp, 'RowNames', opts.geneNames(1:numGenes), 'VariableNames', opts.geneNames));

%======== VISUALIZE RESULTS
close all;
set (gcf, 'Position', [0 0 1200 800], 'Toolbar','None','MenuBar','None');
colormap ('jet');
for g=1:4
    %-------- PLOT EXPERIMENTAL TRAJECTORIES AS HEATMAPS IN (n,t) PLANE --------
    subplot (2, 4, g); hold on;
    imagesc (xntgEXPT(:,:,g)', [0 225]); 
    xlabel ('n'); 
    ylabel ('t'); 
    title ([opts.geneNames{g} ' (experimental)']);
    colorbar ('southoutside');
    %-------- PLOT RE-SIMULATED TRAJECTORIES AS HEATMAPS IN (n,t) PLANE --------
    subplot (2, 4, g+4); hold on;
    imagesc (xntgFIGR(:,:,g)', [0 225]); 
    xlabel ('n'); 
    ylabel ('t'); 
    title ([opts.geneNames{g} ' (resimulated)']);
    colorbar ('southoutside');
    
end
return;






%======== MULTIDIMENSIONAL ARRAY I/O ===============
% In MATLAB R2019a, dlmwrite() has been replaced by writematrix(),
% but there is still no way to write 3D or higher-dimensional arrays
% in human-readable form.
%
% writeArray() and readArray() provide a way to write arbitrary-dimensional
% arrays.
function writeArray (filename, a)
delim = {'\t', '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n'};
%delim = {' ', '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n'};
fmtstr = '%.2f';
dmax = ndims(a);
cmaxd = size(a);
fid = fopen (filename, 'w');
fprintf (fid, '%d dims\n', dmax);
fprintf (fid, '%d elems\n', cmaxd);
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

