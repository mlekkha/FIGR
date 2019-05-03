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
    '=======================================================================\n'...
    '\n']);

%======== Define global structs for options and ODE options 
global opts;
global ODEopts;
global optimopts;
                        

%======== SET OPTIONS (see README.md for description) ========
opts = struct(  'debug', 0, ...                 
                'slopethresh', 1.0, ...         
                'exprthresh', 100.0, ...
                'splinesmoothing', 0.01, ...
                'spatialsmoothing', 0.5, ...
                'minborder_expr_ratio', 0.01, ...
                'Rld_method', 'kink', ...
                'Rld_tsafety', 3, ...
                'synthesisfunction', 'synthesis_sigmoid_sqrt', ...
                'ODEAbsTol', 1e-3, ...
                'ODEsolver', 'ode45');

%======== set the ODE options
ODEopts = odeset('AbsTol', opts.ODEAbsTol);


numGenes = 4;  % 4 genes are "internal"; 3 genes are upstream regulators

%======== Start timer
tic;

%======== READ EXPERIMENTAL TRAJECTORIES xntg(:,1,:) AND TIMEPOINTS tt =======
xntgEXPT = readArray ('fly_xntg.txt');  % xntgEXPT is now a 58x9x7 array
tt       = readArray ('fly_tt.txt');

%======== Print status
disp ('Read data... ');
disp ('Starting FIGR... ');

%======== INFER GRN PARAMETERS grnFIGR
[grnFIGR, diagnostics] = infer (opts, xntgEXPT, tt, numGenes);
yntgEXPT = diagnostics.yntg;

%======== Print status
disp ('FIGR complete... ');
disp ('Starting refinement... ');

%======== REFINE GRN BY NELDER-MEAD ========
xntgFLAT = reshape(xntgEXPT, 522, 7);

% set optimization options for refinement
packed_paramvec = packParams(grnFIGR, numGenes);
optimopts = optimset( 'Display', 'Iter', ...
                      'MaxFunEvals', 200*length(packed_paramvec), ...
                      'MaxIter', 200*length(packed_paramvec));
                        
% check if MEX exists for this platform. Run the mex function if true
% otherwise run interpreted (but slower) MATLAB .m code
if (exist('refineFIGRParams_mex') == 3)
    disp('MEX file for refinement found. Running compiled code...');
    grnREF = refineFIGRParams_mex(grnFIGR, xntgFLAT, tt);
else
    fprintf(1, ['MEX file for refinement ' ...
    '(refineFIGRParams_mex.{mexa64/mexmaci64/mexw64})' ...
    ' not found.\n See README.md for instructions for compiling the ' ...
    ' MEX file.\n\n Running interpreted (but slower) .m code...']);
    grnREF = refineFIGRParams(grnFIGR, xntgFLAT, tt);
end    


%======== Print status
disp ('Nelder-Mead refinement complete... ');

%======== Print time taken
fprintf(1, 'Total time elapsed: %f\n', toc);

%======== RECOMPUTE TRAJECTORIES xntg ========
[xntgREF] = computeTrajs (opts, grnREF, xntgEXPT, tt);

%======== DISPLAY grnSA AND grnFIGR FOR COMPARISON
disp ('GRN parameters inferred using FIGR and refined with Nelder-Mead: ');
disp ('grnREF = ');
disp (struct2table (grnREF));
disp ('grnREF.Tgg = ');
tmp = grnREF.Tgg;

opts.geneNames = ...
            {'Hunchback','Kruppel','Giant','Knirps','Bicoid','Caudal','Tailles'};  % HARDWIRED
disp (array2table (tmp, 'RowNames', opts.geneNames(1:numGenes), 'VariableNames', opts.geneNames));

%======== VISUALIZE RESULTS
close all;
set (gcf, 'Position', [0 0 1200 600], 'Toolbar','None','MenuBar','None');
colormap ('jet');
for g=1:4
    %-------- PLOT EXPERIMENTAL TRAJECTORIES AS HEATMAPS IN (n,t) PLANE --------
    subplot (2, 4, g); hold on;
    imagesc (flipud(xntgEXPT(:,:,g)'), [0 225]); 
    xlabel ('n'); 
    ylabel ('t'); 
    title ([opts.geneNames{g} ' (experimental)']);
    ah = gca;
    set(ah, 'XTick', [10:10:50], 'XTickLabel', [45:10:85]);
    set(ah, 'YTick', [2:2:8], 'YTickLabel', [43.6:-13.5:3.1]);
    if (g == 4)
        oldpos = get(ah, 'Position');
        cbh = colorbar(ah, 'Location', 'manual', ...
          'Position', [oldpos(1)+oldpos(3)+0.01 oldpos(2) 0.01 oldpos(4)]);
    end    
    %-------- PLOT RE-SIMULATED TRAJECTORIES AS HEATMAPS IN (n,t) PLANE --------
    subplot (2, 4, g+4); hold on;
    imagesc (flipud(xntgREF(:,:,g)'), [0 225]); 
    xlabel ('n'); 
    ylabel ('t'); 
    title ([opts.geneNames{g} ' (resimulated)']);
    ah = gca;
    set(ah, 'XTick', [10:10:50], 'XTickLabel', [45:10:85]);
    set(ah, 'YTick', [2:2:8], 'YTickLabel', [43.6:-13.5:3.1]);
    
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

