fprintf ('\n\n');
disp ('&&&&&&&&&&&&&&&&&& main &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp ('RUN THIS MATLAB SCRIPT AS FOLLOWS:              ');
disp ('   opts = defaultOpts;                             ');
disp ('   opts.debug = 2;                                 ');
disp ('   opts.infertype = 1;    % if you wish to change infertype ');
disp ('   main                                                        ');
disp ('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
fprintf ('\n');

% INPUT:    opts struct 
% OUTPUT:   sets many variables in workspace

%======== defaults are set using defaultOpts;
% but individual functions may need to process defaults
% e.g., opts.geneNames='auto' ---->  opts.geneNames={'A','B','C',...}
%======== unfold options include the timepoints file with -j
opts.unfoldopts = [opts.unfoldopts ' -j ' opts.timepoints];

if (opts.debug >= 1) ; disp ('opts = '); disp(opts); end

%======== READ grnSA (Gene Regulatory Network inferred by Simulated Annealing)
% FROM THE FILE hkgn58c14k1_002
[numGenes,numExternals,grnSA,xntgEXPT,nuclei,tt] ...
    = extractData (opts, opts.grnFile);

opts.geneNames = {'H','K','G','N','B','C','T'};  % HARDWIRED


%======== INFER grnCBI (Gene network inferred by Classif-Based Inference)
if (~opts.infertype)   % if infertype==0
    
    %======== grnSA ---> xntgSA
    computeTrajs_fh = str2func(opts.computetrajsfunction);
    [xntgSA] = computeTrajs_fh (opts, grnSA, xntgEXPT, tt);
    
    %======== xntgSA ---> grnCBI
    [grnCBI, diagnostics] = infer (opts, xntgSA, tt, numGenes);
    warning ('YLL: some unsupported stuff ');

else % if infertype==1
    
    %======== xntgEXPT ---> grnCBI
    [grnCBI, diagnostics] = infer (opts, xntgEXPT, tt, numGenes);
    yntgEXPT = diagnostics.yntg;
end

%======== grnCBI ---> xntgCBI
computeTrajs_fh = str2func(opts.computetrajsfunction);
[xntgCBI] = computeTrajs_fh (opts, grnCBI, xntgEXPT, tt);

%======== DISPLAY grnSA AND grnCBI FOR COMPARISON
disp ("grnSA = ");
disp (struct2table (grnSA));
disp ("grnSA.Tgg = ");
tmp = grnSA.Tgg;
disp (array2table (tmp, 'RowNames', opts.geneNames(1:numGenes), 'VariableNames', opts.geneNames));

disp ("grnCBI = ");
disp (struct2table (grnCBI));
disp ("grnCBI.Tgg = ");
tmp = grnCBI.Tgg;
disp (array2table (tmp, 'RowNames', opts.geneNames(1:numGenes), 'VariableNames', opts.geneNames));
%tmp = num2cell (tmp);
%fun = @(x) sprintf('%0.2f', x);
%tmp = cellfun(fun, tmp, 'UniformOutput',0);
%disp (cell2table (tmp, 'RowNames', opts.geneNames(1:numGenes), 'VariableNames', opts.geneNames));
% fprintf(1,'\n    \t R \t Half-life \t Max Expression \t D \n');
% for g=1:numGenes
%     fprintf(1, '%s:\t %f \t %f \t %f \t %f \n', geneNames{g}, ...
%         Rld_inferred(g,1), ...
%         log(2)/Rld_inferred(g,2), ...
%         Rld_inferred(g,1)/Rld_inferred(g,2), ...
%         Rld_inferred(g,3));
% end

%======== COLOR PLOT
disp ('&&&&&&&&&&& NOW MAKE PLOTS USING COMMANDS LIKE:  ');
disp ('   plotGRN (opts, grnSA)  ; title (''grnSA:  T_{gf}, h_g, R_g, \lambda_g, D_g'') ');
disp ('   plotGRN (opts, grnCBI) ; title (''grnCBI: T_{gf}, h_g, R_g, \lambda_g, D_g'') ');

%disp ('   plotGRNComparison (opts,grnSA,grnCBI)           ');
disp ('   plotOnoffstateHeatmap (opts,yntgEXPT,tt,1)           ');
disp ('   plotExpressionHeatmap (opts,xntgEXPT,tt,1)                      ');
disp ('   for g=1:7; plotExpressionHeatmap (opts,xntgEXPT,tt,g);pause(0.5);end   ');
disp ('   plot2DProjections (opts,xntgEXPT,yntgEXPT,tt,1,3,6)       ');
disp ('   for n=1:58; plotTrajs (opts,xntgEXPT,xntgSA,tt,n); pause(0.5); end  '); %nucleus 50
disp ('   plotExpression3D(opts,xntgEXPT,tt,yntgEXPT,1)            ');
disp ('   plotExpression3DComparison(opts,xntgEXPT,xntgCBI,tt,2)           ');
disp ('             ');
disp (' WE ARE INTERESTED IN COMPARING xntgCBI WITH xntgEXPT            ');
disp ('   plotTrajs (opts,xntgEXPT,xntgCBI,tt,3)          ');
disp ('   plotExpressionHeatmap (opts,xntgEXPT,tt,1) ; title ''EXPT''     ');
disp ('   plotExpressionHeatmap (opts,xntgCBI,tt,1)  ; title ''CBI''        ');
disp ('             ');
disp (' NOT YET PROPERLY IMPLEMENTED:  ');
disp ('   plotTrajsComparison(...)  -- should also calc andprint chi^2    ');
disp ('             ');
disp (' YOU MAY WISH TO DO close all OR showfigs AT SOME STAGE ');
disp (' YOU MAY WISH TO PREPARE A SCRIPT SUCH AS plt.m    ');
disp ('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');



clear j param_names defaultParams computeTrajs_fh tmp
return;





% See runToyModel.m and runToyModels.m for other top-level functions.
