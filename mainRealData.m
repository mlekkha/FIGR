%======== This is a Matlab script! ========

fprintf ('\n\n\n\n');
fprintf ('=======================================\n');

%======== READ FLY DATA FILE hkgn58c14k1_002 INTO Rg,lambdag,etc.
% Call extractData() but prevent it from spitting out its debugging information
evalc ('[numGenes,numExternals,Rg,lambdag,Tgg,Xntg,nuclei,tt,geneNames] = extractData (''~/genecircuits/hkgn58c14k1_002'') ');

%======== COMPUTE TRAJECTORIES xntg BY CALLING unfold
[xntg] = computeTrajsUsingUnfold (Rg,lambdag,Tgg,Xntg,nuclei,tt,geneNames);

%======== PLOT TRAJECTORIES OF EACH NUCLEUS, OR EVOLUTION OF SPATIAL PROFILES
%plotTrajs (Xntg,xntg,nuclei,tt,geneNames);
%plotSpatial (Xntg,xntg,nuclei,tt,geneNames);

%======== INFER TggInferred - ----------  FROM REAL DATA!!!!!!!!!!!
%                     ADMITTEDLY, DOESN'T WORK AS WELL,
%                BUT MAYBE BECAUSE I'M COMPARING WITH T-MATRIX PARS OPTIMIZED FOR SIGMOID
%                WHICH DON'T REALLY FIT THE REAL DATA WITH HEAVISIDE MODEL???

%TggInferred = infer (numGenes,numExternals,xntg,nuclei,tt);
TggInferred = infer (numGenes,numExternals,Xntg,nuclei,tt);




%======== COMPUTE TRAJECTORIES xntg BY CALLING unfold
%%%% THIS IS NEW. 2018-9-25.
[xntgFromCBIParams] = computeTrajsUsingUnfold (Rg,lambdag,TggInferred,Xntg,nuclei,tt,geneNames);
plotSpatial (Xntg,xntgFromCBIParams,nuclei,tt,geneNames);


%======== NORMALIZE AND COMPARE
for g=1:numGenes
    Tgg(g,:) = Tgg(g,:) / norm(Tgg(g,:));
    TggInferred(g,:) = TggInferred(g,:) / norm(TggInferred(g,:));
end
Tgg = Tgg * 1000;
TggInferred = TggInferred * 1000;

%======== COLOR PLOT
figure (3);
myColormap = flipud (redbluecmap (21));
myColormap = myColormap (:, [1 3 2]);  % leave R but swap G and B components
colormap (myColormap);

subplot (2,1,1);
imagesc (Tgg, [-100 100]);
colorbar;
%axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:8],'xticklabel', {'H','K','G','N','B','C','T','h'});
set(gca,'ytick',[1:4],'yticklabel', {'H','K','G','N'});
title ('Theoretical T-matrix and h-vector');
for g=1:size(Tgg,1)
    for f=1:size(Tgg,2)
        text (f, g, num2str(round(Tgg(g,f) )));
    end
end

subplot (2,1,2);
imagesc (TggInferred, [-100 100]);
colorbar;
%axis image; % force 1:1 aspect ratio
set(gca,'xtick',[1:8],'xticklabel', {'H','K','G','N','B','C','D','h'});
set(gca,'ytick',[1:4],'yticklabel', {'H','K','G','N'});
title ('Inferred T-matrix and h-vector');
for g=1:size(Tgg,1)
    for f=1:size(Tgg,2)
        text (f, g, num2str(round(TggInferred(g,f) )));
    end
end
print ('fig3', '-dpng');

%======== T AND TInferred VERSUS 32 GENE INDICES
figure (4); clf; hold on;
title ( sprintf ('T-matrix elements') );
xlabel ('Element no.');
ylabel ('Element');
xlim ([0 numel(Tgg)])
ylim ([-100 200]);
grid on;
plot (Tgg(:), 'bx', 'LineWidth', 3);
plot (TggInferred(:), 'ro', 'LineWidth', 3);
legend ('Actual', 'Inferred');
x = linspace (-1000,1000); % Manually draw axes
y = linspace (0,0);
plot(x,y,'k-');
print ('fig4', '-dpng');

%======== COVARIANCE PLOT
figure (5); clf; hold on;
title ( sprintf ('Covariance of T vs Tinferred') );
xlabel ('T_actual');
ylabel ('T_inferred');
xlim ([-100 200]);
ylim ([-100 200]);
grid on;
plot (Tgg(:), TggInferred(:), 'ro', 'LineWidth', 3);
x = linspace (-1000,1000); % Manually draw axes
y = linspace (0,0);
plot(x,y,'k-');
plot(y,x,'k-');
print ('fig5', '-dpng');

dlmwrite ('~/genecircuits/TggManu.dat', Tgg, 'delimiter', '\t', 'precision', 2);
dlmwrite ('~/genecircuits/TggInferred.dat', TggInferred, 'delimiter', '\t', 'precision', 2);

return;
