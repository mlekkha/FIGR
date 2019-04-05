function [opts] = defaultOpts ()

%======== DEFINE default values of parameters
opts.debug = 0;
opts.circdir = '~/genecircuits/';  %CURRENTLY NOT BEING USED PROPERLY
opts.datfile = '~/genecircuits/Xntg.dat';
opts.unfold = '/usr/local/fly/fly/unfold';  %HARDWIRED THE FOLDER NAME
opts.timepoints = '~/genecircuits/timepoints';
opts.tempcirc = '~/genecircuits/networkPars';
opts.unfoldout = '~/genecircuits/trajectories';
opts.grnFile = '~/genecircuits/hkgn58c14k1_002';
opts.TggTheorFile = 'TggManu.dat';
opts.TggInferFile = 'TggInferred.dat';

% parameters for the inference
opts.slopethresh = 1.0;  % YL 2018-11-6
opts.exprthresh = 100.0; % YL 2018-11-6
opts.splinesmoothing = 0.01;
opts.spatialsmoothing = 0.5;
opts.minborder_expr_ratio = 0.01;
opts.Rld_method = 'slope';
opts.Rld_tsafety = 3;
opts.Rld_thresh_on_frac = 0.2;
opts.Rld_thresh_off_frac = 0.01;

% parameters for computeTrajs
opts.unfoldopts = '-s r4 -i 0.01 -a 0.001 -x input -g h ';
opts.synthesisfunction = 'synthesis_sigmoid_sqrt';
opts.computetrajsfunction = 'computeTrajs';

% ODE solver parameters
opts.ODEAbsTol = 1e-2;
opts.ODEsolver = 'ode45';

% inference type: 0: simulated data, 1: read data
opts.infertype = 0;
opts.unfoldopts = [opts.unfoldopts ' -j ' opts.timepoints];

%opts.geneNames = geneNames;
end
