%============================================================
% Modified Principal Component Analysis for Binary-Classified Data
% Yen Lee Loh 2020-6-7
%
% Demonstrates the use of mpca().
%============================================================
function [] = mpcaDemo ()
clc; fprintf ('\n============== mpcaDemo =====================\n');

fprintf ('1. Toy model with datapoints in elliptical disk \n');
fprintf ('2. Toy model with datapoints in ellipsoid \n');
fprintf ('3. Fly dataset (HKGNBCT) \n');
fprintf ('4. ERYNEU dataset (12 genes + 1 external regulator) \n');
choice = input ("Which of the above examples do you want to run (Enter a number)? ");
if (choice==1)
    %======== EXAMPLE 1 ==================================================
    %======== PREPARE DATAPOINTS xkg AND CLASSES yk: SYNTHETIC GAUSSIAN DATA
    fprintf ("=============== EXAMPLE 1 README!!!!!!!!!!! ======= \n");
    fprintf ("This example shows data points where gene G is classified as ON (green) or OFF (red)\n");
    fprintf ("  depending on expressions of genes A, B, C.\n");
    fprintf ("The data points form a disk,\n");
    fprintf ("  and they are cleanly separated in the projections.\n");
    fprintf ("  (For more general data the projections will overlap.)\n");
    fprintf ("Here gene C activates gene G, whereas A and B have almost no effect.\n");
    fprintf ("\n");
    fprintf ("The MPCA algorithm chooses principal direction 1 parallel to \n");
    fprintf ("  the decision boundary normal, T.  This is close to the C direction.\n");
    fprintf ("MPCA chooses principal direction 2 close to the B direction\n");
    fprintf ("  so that the points are spread out as much as possible.\n");
    fprintf ("Principal direction 3 is the normal to the disk\n");
    fprintf ("  in this example.\n");
    
    rng ('default');
    kmax = 400;                  % number of data points
    gmax = 3;                    % dimensionality
    geneNames = ["A" "B" "C"]';
    center = [8.  5.  7.];       % center of multidim Gaussian distribution
    widths = [4.  3.  .1];       % widths of multidim Gaussian distribution
    widthOfSigmoid = 0.2;        % width of sigmoid function
    offset = -3.0;               % imbalance parameter
    rotationMatrix = rot3d (-90*pi/180.0, [0 1 0]) ...
        *  rot3d (15*pi/180.0, [1 1 1]);
    % rotationMatrix = randOrthMat(gmax);
    xkg = NaN (kmax, gmax);
    yk  = NaN (kmax, 1);
    for k=1:kmax
        while (true)   % restricdt to nice flat ellipse
            xkg(k,:) = normrnd (0, widths, [1 gmax]);
            if (norm(xkg(k,:) ./ widths) < 1.5); break; end
        end
        distFromBoundary = xkg(k,1) - offset;
        probOfBeingOn = 1 / (1 + exp(-distFromBoundary/widthOfSigmoid));
        yk(k) = sign (probOfBeingOn - rand());
    end
    xgk = xkg';
    xgk = rotationMatrix * xgk;  % left-multiply active rotation onto dataset
    xkg = xgk';
    xkg = xkg + center;          % shift
    
    imageSize = 14.;
    scal = imageSize*.8;
    
elseif (choice==2)
    
    %======== EXAMPLE 2 ==================================================
    %======== PREPARE DATAPOINTS xkg AND CLASSES yk: SYNTHETIC GAUSSIAN DATA
    fprintf ("=============== EXAMPLE 2 README!!!!!!!!!!! ======= \n");
    fprintf ("This example shows data points where gene G is classified as ON (green) or OFF (red)\n");
    fprintf ("  depending on expressions of genes A, B, C.\n");
    fprintf ("The data points form an ellipsoid.\n");
    fprintf ("  They are NOT cleanly separated in the projections.\n");
    fprintf ("\n");
    fprintf ("Nevertheless, the MPCA plot illustrates that the data are cleanly separated\n");
    fprintf ("  by the decision boundary (click and drag to rotate to verify this)!\n");
    fprintf ("\n");
    fprintf ("Rotate from the 1-2 plane to the 1-3 plane to verify that\n");
    fprintf ("  the data have greater spread in the 2 direction.\n");
    fprintf ("\n");
    
    rng ('default');
    kmax = 400;                  % number of data points
    gmax = 3;                    % dimensionality
    geneNames = ["A" "B" "C"]';
    center = [8.  5.  7.];       % center of multidim Gaussian distribution
    widths = [4.  3.  8.];       % widths of multidim Gaussian distribution
    widthOfSigmoid = 0.2;        % width of sigmoid function
    offset = -3.0;               % imbalance parameter
    rotationMatrix = rot3d (80*pi/180.0, [1 1 1]);
    xkg = NaN (kmax, gmax);
    yk  = NaN (kmax, 1);
    for k=1:kmax
        while (true)
            xkg(k,:) = normrnd (0, widths, [1 gmax]);
            if (norm(xkg(k,:) ./ widths) < 1.5); break; end  % make a nice ellipsoid
        end
        distFromBoundary = xkg(k,1) - offset;
        probOfBeingOn = 1 / (1 + exp(-distFromBoundary/widthOfSigmoid));
        yk(k) = sign (probOfBeingOn - rand());
    end
    xgk = xkg';
    xgk = rotationMatrix * xgk;  % left-multiply active rotation onto dataset
    xkg = xgk';
    xkg = xkg + center;          % shift
    
    imageSize = 14.;
    scal = imageSize*.8;
    
    
elseif (choice==3)
    %======== EXAMPLE 3 ==================================================
    % PREPARE DATAPOINTS xkg AND CLASSES yk: FLY DATA    
    geneNames = ["H" "K" "G" "N" "B" "C" "T"]';
    xntg = load('xntg.mat').xntg;       % note syntax
    yntg = load('yntg.mat').yntg;
    [nmax tmax gmax] = size (xntg);     % number of nuclei, timepts, genes
    kmax = nmax*tmax;                   % number of datapoints
    xkg = reshape (xntg, [kmax gmax]);
    ykg = reshape (yntg, [kmax gmax]);
    
    fprintf ("The order of genes is:\n"); disp ( [(1:numel(geneNames))' geneNames] ); 
    gTarget = input ("Which gene's ON/OFF states to plot (enter a number)? ");
    
    yk = ykg(:,gTarget);
    
    imageSize = 250.;
    scal = imageSize*.8;
    
elseif (choice==4)
    %======== EXAMPLE 4 ==================================================
    % PREPARE DATAPOINTS xkg AND CLASSES yk: ESPER DATA

    [xntg tt nucleusNames geneNames] = loadNCBIFiles("xntg.txt", "tn.txt");
    [yntg] = loadMDA ("eryneu_yntg.mda");
    [nmax tmax gmax] = size (xntg);     % number of nuclei, timepts, genes
    kmax = nmax*tmax;                   % number of datapoints
    xkg = reshape (xntg, [kmax gmax]);
    ykg = reshape (yntg, [kmax gmax]);
    
    fprintf ("The order of genes is:\n"); disp ( [(1:numel(geneNames))' geneNames] ); 
    gTarget = input ("Which gene's ON/OFF states to plot (enter a number)? ");
    
    yk = ykg(:,gTarget);
    
    imageSize = 1.;
    scal = imageSize*.8;
   
else
    fprintf ("Invalid choice!\n");
    return;
end



%======== FIND CLASSIFICATION HYPERPLANE USING LOGISTIC REGRESSION
% GIVEN xkg AND yk, FIND Tg AND h
beta = glmfit (xkg, max(yk,0), 'binomial', 'link', 'logit');  % We get h,T1,T2,...,
Tg   = beta(2:end);
h    = beta(1);

fprintf ("Number of datapoints kmax = %d   (%d ON, %d OFF) \n", kmax, sum(yk>0), sum(yk<0));
fprintf ("Dimensionality       gmax = %d\n", gmax);
fprintf ('Classification boundary normal Tg = \n'); disp (Tg);

%======== PERFORM MODIFIED PRINCIPAL COMPONENT ANALYSIS (MPCA)
% Columns of UggPrime are principal directions
[UggPrime, xkgPrime] = mpca (xkg, yk, Tg);
UgPrimeg = inv(UggPrime); % transformation from primed coords back to original coords







%============================================================
% VISUALIZATION SECTION
%============================================================
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontWeight','Normal');  % Bold
figure(1); clf; set (gcf, 'Position', [10 10 1200 900]);
colorON  = [0 .6 0];  % green
colorOFF = [.8 0 0];  % red
colors = NaN (kmax, 3);
for k=1:kmax
    if yk(k)>0 ; colors(k,:) = colorON;
    else       ; colors(k,:) = colorOFF;
    end
end

handleR3D = rotate3d; handleR3D.Enable = 'on'; rotate3d on;
%======== VISUALIZE ORIGINAL DATA AS 3D SCATTER PLOT
subplot (2,3,1); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter3 (xkg(:,1), xkg(:,2), xkg(:,3), 20, colors);
title (sprintf ('3D Plot (Rotate Me!)'));
view (30, 15);
arrow3 (scal*[1 0 0], geneNames(1), 'black'); xlabel (geneNames(1));
arrow3 (scal*[0 1 0], geneNames(2), 'black'); ylabel (geneNames(2));
arrow3 (scal*[0 0 1], geneNames(3), 'black'); zlabel (geneNames(3));
xlim ([-imageSize imageSize]);
ylim ([-imageSize imageSize]);
zlim ([-imageSize imageSize]);
%xticks (-15:5:15); yticks (-15:5:15); zticks (-15:5:15);

%======== VISUALIZE DATA IN "PRINCIPAL" FRAME (3D SCATTER PLOT)
subplot (2,3,3); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter3 (xkgPrime(:,1), xkgPrime(:,2), xkgPrime(:,3), 20, colors);
title (sprintf ('Principal Coordinates (Rotate Me!)'));
view (30, 15);
arrow3 (scal*[1 0 0], "1", 'red'); xlabel ('PD1');
arrow3 (scal*[0 1 0], "2", 'red'); ylabel ('PD2'); 
arrow3 (scal*[0 0 1], "3", 'red'); zlabel ('PD3');
for f=1:gmax  % loop over transcription factors f that INFLUENCE gene g
    arrow3 (scal*UgPrimeg([1 2 3], f), geneNames(f), 'black');
end
xlim ([-imageSize imageSize]);
ylim ([-imageSize imageSize]);
zlim ([-imageSize imageSize]);
%---- Plot 2d decision boundary (line along y' axis)
%line ([0 0], [-imageSize imageSize], 'LineStyle', '--', 'Color', 'black');
%---- Plot 3d decision boundary (y'-z' plane)
x0 = -h/norm(Tg);               % distance of decision boundary from origin
patch (...
    [x0 x0 x0 x0], ...          % x' coordinates of rectangle
    imageSize*[1 -1 -1 1], ...  % y' coordinates of rectangle
    imageSize*[1 1 -1 -1], ...  % z' coordinates of rectangle
    [.5 .5 .5], 'FaceAlpha', .3);


%======== VISUALIZE ORIGINAL DATA AS 2D SCATTER PLOT (A-B PLANE etc.)
g1 = 1; g2 = 2;
subplot (2,3,4); 
cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
scatter (xkg(:,g1), xkg(:,g2), 20, colors);
arrow (scal*[1 0], geneNames(g1), 'black');
arrow (scal*[0 1], geneNames(g2), 'black');
arrow (scal*UggPrime([1 2],g1), "1", 'red');
arrow (scal*UggPrime([1 2],g2), "2", 'red');
title ("2D Projection");
xlabel (geneNames(g1)); ylabel (geneNames(g2));
setAllowAxesRotate (handleR3D, gca, false);  % don't allow rotating THIS

g1 = 1; g2 = 3;
subplot (2,3,5); 
cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
scatter (xkg(:,g1), xkg(:,g2), 20, colors);
arrow (scal*[1 0], geneNames(g1), 'black');
arrow (scal*[0 1], geneNames(g2), 'black');
arrow (scal*UggPrime([1 2],g1), "1", 'red');
arrow (scal*UggPrime([1 2],g2), "2", 'red');
title ("2D Projection");
xlabel (geneNames(g1)); ylabel (geneNames(g2));
setAllowAxesRotate (handleR3D, gca, false);  % don't allow rotating THIS

g1 = 2; g2 = 3;
subplot (2,3,6); 
cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
scatter (xkg(:,g1), xkg(:,g2), 20, colors);
arrow (scal*[1 0], geneNames(g1), 'black');
arrow (scal*[0 1], geneNames(g2), 'black');
arrow (scal*UggPrime([1 2],g1), "1", 'red');
arrow (scal*UggPrime([1 2],g2), "2", 'red');
title ("2D Projection");
xlabel (geneNames(g1)); ylabel (geneNames(g2));
setAllowAxesRotate (handleR3D, gca, false);  % don't allow rotating THIS
return;

%============================================================
% NESTED UTILITY FUNCTIONS BELOW
%============================================================
    function arrow (dir, label, color)
        quiver (0, 0, dir(1), dir(2), color, 'LineWidth', 1.5);
        text (1.05*dir(1), 1.05*dir(2), label, 'FontSize', 16, 'FontWeight', 'Bold');
    end
    function arrow3 (dir, label, color)
        quiver3 (0, 0, 0, dir(1), dir(2), dir(3), color, 'LineWidth', 1.5);
        text (1.05*dir(1), 1.05*dir(2), 1.05*dir(3), label, 'FontSize', 16, 'FontWeight', 'Bold');
    end

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
