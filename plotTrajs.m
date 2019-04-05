function [] = plotTrajs (opts,Xntg,xntg,tt,n) % n is nucleus index

fprintf ('================ plotTrajs() =======================\n');
disp ('THIS FUNCTION IS A BIT AMBIGUOUS');
disp ('IN THAT IT PLOTS EXPT TRAJ AS WELL AS SIMULATED TRAJ....');
disp ('for certain values infertype   this doesn not make sense');

numNuclei     = size (Xntg, 1);
numTimepoints = size (Xntg, 2);
numRegulators = size (Xntg, 3);
%numExternals  = numRegulators - numGenes;
geneNames = opts.geneNames;

numGenes = numRegulators;  % JUST PLOT EVERYTHING

%======= Resistor color codes for 1234567890
colorList = .2 * [3 2 1; 5 0 0; 5 3 0; 5 5 0; 0 5 0; 0 0 5; 5 0 5; 4 4 4; 0 5 5; 0 0 0];



figure (2); clf; hold on;
title ( sprintf ('Trajectories x(n=%d,t,g)', n ), 'FontSize', 14 );
xlabel ('Time t');
ylabel ('Concentrations x');
%ylim ([0 250]);
%-------- PLOT EXPT TRAJECTOREIS OF HKGN
for g=1:numGenes
    theColor = colorList(g,:);  % '--x'
    plot (tt, squeeze(Xntg(n,:,g)), 'x', ...
        'DisplayName', [geneNames{g} ' expt'], ...
        'LineWidth', 1, 'Color', theColor, 'MarkerFaceColor', theColor);
end
%-------- PLOT SIMULATED TRAJECTORIES OF HKGN
for g=1:numGenes
    theColor = colorList(g,:);
    plot (tt, squeeze(xntg(n,:,g)), '-', ...
        'DisplayName', [geneNames{g} ' sim'], ...
        'LineWidth', 3, 'Color', theColor, 'MarkerFaceColor', theColor);
end
legend ('Position',[0.9 0.5 0.1 0.2]); % (geneNames{1:numGenes});
%-------- PLOT SPLINE FIT TO SIMULATED TRAJS
for g=1:numGenes
    theColor = colorList(g,:);
    tFineGrid = linspace (min(tt), max(tt), 501);    
    mySpline = csaps (tt, xntg(n,:,g), opts.splinesmoothing);   % WE CHECKED THAT WE NEED THIS MUCH SMOOTHING
    mySplDer = fnder (mySpline);         % Calculate derivative w.r.t. gene g!
    plot (tFineGrid, fnval (mySpline, tFineGrid), 'Color', theColor, ...
                'DisplayName', [geneNames{g} ' smoothed']);
end




% figure (3); clf; hold on;
% %-------- COMPARE SLOPES
% for g=1:numGenes
%     theColor = colorList(g,:);
%     tFineGrid = linspace (min(tt), max(tt), 501);
%     
%     mySpline = csaps (tt, xntg(n,:,g), opts.splinesmoothing);   % WE CHECKED THAT WE NEED THIS MUCH SMOOTHING
%     mySplDer = fnder (mySpline);         % Calculate derivative w.r.t. gene g!
%     plot (tFineGrid, fnval (mySplDer, tFineGrid), '-.', 'Color', theColor);
% end 


fprintf ('=======================================\n\n\n\n');
end

