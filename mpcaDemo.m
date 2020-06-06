%============================================================
% Modified Principal Component Analysis for Binary-Classified Data
% Yen Lee Loh 2020-6-5
% 
% Demonstrates the use of mpca().
% THE WHOLE THING IS MESSY BECAUSE MATLAB LIKES COLUMN VECTORS
% AND x(k,g) CONSISTS OF ROW VECTORS
%============================================================
function [] = mpcaDemo () 
clc; fprintf ('\n============== mpcaDemo =====================\n');

%======== PREPARE DATAPOINTS xkg AND CLASSES yk: FLY DATA
% (NOT READY TO UNCOMMENT YET --- TOO MANY DIMENSIONS)
%
% geneNames = {'H', 'K', 'G', 'N', 'B', 'C', 'T'}; 
% xntg = load('xntg.mat').xntg;       % note syntax
% yntg = load('yntg.mat').yntg;
% [nmax tmax gmax] = size (xntg);     % number of nuclei, timepts, genes
% kmax = nmax*tmax;                   % number of datapoints
% xkg = reshape (xntg, [kmax gmax]);
% ykg = reshape (yntg, [kmax gmax]);
% gTarget = 2;                % consider ON/OFF state of gene gTarget
% yk = ykg(:,gTarget);

%======== PREPARE DATAPOINTS xkg AND CLASSES yk: SYNTHETIC GAUSSIAN DATA
rng ('default');
kmax = 400;                  % number of data points
gmax = 3;                    % dimensionality
geneNames = ["A" "B" "C"]'; 
center = [8.  5.  7.];       % center of multidim Gaussian distribution
widths = [4.  1.  3.];       % widths of multidim Gaussian distribution
widthOfSigmoid = 0.2;        % width of sigmoid function
offset = -3.0;               % imbalance parameter
%rotationMatrix = rot3d (0*pi/180.0, [0 0 1]); % rotate the multidim Gaussian distribution
rotationMatrix = rot3d (95*pi/180.0, [0 9 1]); % rotate the multidim Gaussian distribution
% rotationMatrix = randOrthMat(gmax);

xkg = NaN (kmax, gmax);
yk  = NaN (kmax, 1);
for k=1:kmax
    xkg(k,:) = normrnd (0, widths, [1 gmax]);
    distFromBoundary = xkg(k,1) - offset;
    probOfBeingOn = 1 / (1 + exp(-distFromBoundary/widthOfSigmoid));
    yk(k) = sign (probOfBeingOn - rand());
end
xkg = xkg * rotationMatrix' + center;  % note transpose

fprintf ("Number of datapoints kmax = %d   (%d ON, %d OFF) \n", kmax, sum(yk>0), sum(yk<0));
fprintf ("Dimensionality       gmax = %d\n", gmax);

%======== FIND CLASSIFICATION HYPERPLANE USING LOGISTIC REGRESSION
beta = glmfit (xkg, max(yk,0), 'binomial', 'link', 'logit');  % We get h,T1,T2,...,
Tg   = beta(2:end);
h    = beta(1);
fprintf ('Classification boundary normal Tg = \n'); disp (Tg);

%======== PERFORM MODIFIED PRINCIPAL COMPONENT ANALYSIS (MPCA)
% basis contains ROW vectors that are the principal directions!
[basis, xkgRot] = mpca (xkg, yk, Tg, h);

fprintf ("==================== READ THIS EXPLANATION! ======= \n");
fprintf ("This example has gene G classified ON (green) and OFF (red) data points.\n");
fprintf ("The data points are NOT cleanly separated in the projections.\n");
fprintf ("However, we see that increasing C tends to switch G OFF; \n");
fprintf ("  increasing B has a slight tendency to switch G ON; \n");
fprintf ("  and A hardly affects the state of G. \n");
fprintf ("\n");
fprintf ("The 'MPCA' plot shows the data rotated to show a clean separation.\n");
fprintf ("The hyperplane is vertical; its normal is horizontal.\n");
fprintf ("The 'MPCA' plot shows that C is almost antiparallel to the hyperplane normal, \n");
fprintf ("  B only has a small projection in the direction of that normal,  \n");
fprintf ("  and A is almost parallel to the hyperplane.  \n");
fprintf ("Also, the 'MPCA' plot chooses direction 2 to be the direction of\n");
fprintf ("  largest spread [with the constraint that 2 is perp to 1]. \n");



%============================================================
% VISUALIZATION SECTION
%============================================================
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontWeight','Normal');  % Bold
figure(1); clf; set (gcf, 'Position', [10 10 1200 900]);
imageSize = 16.;
colorON  = [0 .6 0];  % green
colorOFF = [.8 0 0];  % red
colors = NaN (kmax, 3);
for k=1:kmax
    if yk(k)>0 ; colors(k,:) = colorON;
    else       ; colors(k,:) = colorOFF;
    end
end

%======== VISUALIZE ORIGINAL DATA AS 3D SCATTER PLOT
subplot (2,3,1); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
ind = find (yk>0);
plot3 (xkg(ind,1),xkg(ind,2),xkg(ind,3), 'o', 'Color', colorON);
ind = find (yk<0);
plot3 (xkg(ind,1),xkg(ind,2),xkg(ind,3), 'o', 'Color', colorOFF);
title (sprintf ('3D Plot (Rotate Me!)'));
xlabel ('Gene A'); ylabel ('Gene B'); zlabel ('Gene C');
xlim ([-imageSize imageSize]);
ylim ([-imageSize imageSize]);
zlim ([-imageSize imageSize]);

%======== VISUALIZE ORIGINAL DATA AS 2D SCATTER PLOT (A-B PLANE etc.)
subplot (2,3,4); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter (xkg(:,1), xkg(:,2), 20, colors);
xlabel ('Gene A'); ylabel ('Gene B');
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
title (sprintf ('Projection onto AB Plane'));

subplot (2,3,5); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter (xkg(:,1), xkg(:,3), 20, colors);
xlabel ('Gene A'); ylabel ('Gene C');
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
title (sprintf ('Projection onto AC Plane'));

subplot (2,3,6); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter (xkg(:,2), xkg(:,3), 20, colors);
xlabel ('Gene B'); ylabel ('Gene C');
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
title (sprintf ('Projection onto BC Plane'));

%======== VISUALIZE DATA IN "PRINCIPAL" FRAME
subplot (2,3,3); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter (xkgRot(:,1), xkgRot(:,2), 20, colors);
xlabel ('Direction 1 (Hyperplane Normal T)');
ylabel ('Direction 2 (Largest Spread)');
title (sprintf ('Principal Coordinates'));
for f=1:gmax  % loop over transcription factors f that INFLUENCE gene g
    quiver (0, 0, imageSize*basis(f,1), imageSize*basis(f,2), 'LineWidth', 1.5);
    text (1.05*imageSize*basis(f,1), 1.05*imageSize*basis(f,2), sprintf ('%s', geneNames{f}), 'FontSize', 16, 'FontWeight', 'Bold');
end
xlim ([-imageSize imageSize]);
ylim ([-imageSize imageSize]);
line ([0 0], [-imageSize imageSize], 'LineStyle', '--', 'Color', 'black');

%======== COMPARE WITH STANDARD PCA?
% figure(3); clf; hold on;
% [coeff,score,latent] = pca (xkg);
% biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'H','K','G','N','B','C','T'});

    function R = rot3d (rotAngle, rotAxis)
        s = sin(rotAngle);
        c = cos(rotAngle);
        u = rotAxis(:);
        u = u ./ sqrt(u.' * u);        
        x  = u(1);
        y  = u(2);
        z  = u(3);
        mc = 1 - c;
        R  = [c + x * x * mc,    x * y * mc - z * s,   x * z * mc + y * s; ...
            x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
            x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];        
    end

    function M = randOrthMat(n, tol)
        %=================================================
        % M = RANDORTHMAT(n)
        % generates a random n x n orthogonal real matrix.
        %
        % M = RANDORTHMAT(n,tol)
        % explicitly specifies a thresh value that measures linear dependence
        % of a newly formed column with the existing columns. Defaults to 1e-6.
        %
        % In this version the generated matrix distribution *is* uniform over the manifold
        % O(n) w.r.t. the induced R^(n^2) Lebesgue measure, at a slight computational
        % overhead (randn + normalization, as opposed to rand ).
        %
        % (c) Ofek Shilon , 2006.
        if nargin==1
            tol=1e-6;
        end
        
        M = zeros(n); % prealloc
        
        % gram-schmidt on random column vectors
        
        vi = randn(n,1);
        % the n-dimensional normal distribution has spherical symmetry, which implies
        % that after normalization the drawn vectors would be uniformly distributed on the
        % n-dimensional unit sphere.
        M(:,1) = vi ./ norm(vi);
        
        for i=2:n
            nrm = 0;
            while nrm<tol
                vi = randn(n,1);
                vi = vi -  M(:,1:i-1)  * ( M(:,1:i-1).' * vi )  ;
                nrm = norm(vi);
            end
            M(:,i) = vi ./ nrm;
        end %i
        
    end  % RandOrthMat

end
