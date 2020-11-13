%============================================================
% Principal Component Analysis
% Yen Lee Loh 2020-6-6
% 
% Demonstrates the use of pca().
%============================================================
function [] = pcaDemo () 
clc; fprintf ('\n============== mpcaDemo =====================\n');

%======== PREPARE DATAPOINTS xkg AND CLASSES yk: SYNTHETIC GAUSSIAN DATA
rng ('default');
kmax = 900;                  % number of data points
gmax = 3;                    % dimensionality
geneNames = ["A" "B" "C"]'; 
center    = [4.  4.  4.];    % center of multidim Gaussian distribution
widths    = [6.  3.  0.25];   % widths of multidim Gaussian distribution
%rotationMatrix = rot3d (30*pi/180.0, [0 0 1]); % rotate the multidim Gaussian distribution
rotationMatrix = rot3d (60*pi/180.0, [1 1 1]); % rotate the multidim Gaussian distribution

xkg = NaN (kmax, gmax);
yk  = NaN (kmax, 1);
for k=1:kmax
    xkg(k,:) = normrnd (0, widths, [1 gmax]);
end
xgk = xkg';
xgk = rotationMatrix * xgk;  % left-multiply active rotation onto dataset
xkg = xgk';
xkg = xkg + center;          % shift

fprintf ("Number of datapoints kmax = %d\n", kmax);
fprintf ("Dimensionality       gmax = %d\n", gmax);

%======== PERFORM PRINCIPAL COMPONENT ANALYSIS (PCA)
% Columns of princDirs are principal directions
[UggPrime, xkgPrime] = pca (xkg);
UgPrimeg = inv(UggPrime); % transformation from primed coords back to original coords

%============================================================
% VISUALIZATION SECTION
%============================================================
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontWeight','Normal');  % Bold
figure(1); clf; set (gcf, 'Position', [10 10 1200 650]);
imageSize = 17.;
scal = imageSize*.8;
myColor = [.8 0 0];  % red
colors = repmat (myColor, [kmax 1]);

%======== VISUALIZE ORIGINAL DATA AS 3D SCATTER PLOT
subplot (2,3,1); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
plot3 (xkg(:,1),xkg(:,2),xkg(:,3), 'o', 'Color', myColor);
title (sprintf ('3D Plot (Rotate Me!)'));
xlabel ('Gene A'); ylabel ('Gene B'); zlabel ('Gene C');
xlim ([-imageSize imageSize]);
ylim ([-imageSize imageSize]);
zlim ([-imageSize imageSize]);
view (30, 15);


%======== VISUALIZE ORIGINAL DATA AS 2D SCATTER PLOT (A-B PLANE etc.)
subplot (2,3,4); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
xticks (-15:5:15); yticks (-15:5:15); 
scatter (xkg(:,1), xkg(:,2), 20, colors);
arrow (scal*[1 0], "A");
arrow (scal*[0 1], "B");
arrow (scal*UggPrime([1 2],1), "1");
arrow (scal*UggPrime([1 2],2), "2");
title ("Projection onto AB Plane");
xlabel ("Gene A"); ylabel ("Gene B");

subplot (2,3,5); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
xticks (-15:5:15); yticks (-15:5:15); 
scatter (xkg(:,1), xkg(:,3), 20, colors);
arrow (scal*[1 0], "A");
arrow (scal*[0 1], "C");
arrow (scal*UggPrime([1 3],1), "1");
arrow (scal*UggPrime([1 3],2), "2");
title ("Projection onto AC Plane");
xlabel ("Gene A"); ylabel ("Gene C");

subplot (2,3,6); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
xlim ([-imageSize imageSize]); ylim ([-imageSize imageSize]);
xticks (-15:5:15); yticks (-15:5:15); 
scatter (xkg(:,2), xkg(:,3), 20, colors);
arrow (scal*[1 0], "B");
arrow (scal*[0 1], "C");
arrow (scal*UggPrime([2 3],1), "1");
arrow (scal*UggPrime([2 3],2), "2");
title ("Projection onto BC Plane");
xlabel ("Gene B"); ylabel ("Gene C");


%======== VISUALIZE DATA IN "PRINCIPAL" FRAME
subplot (2,3,3); cla; hold on; grid on; set (gca,'DataAspectRatio', [1 1 1]);
scatter (xkgPrime(:,1), xkgPrime(:,2), 20, colors);
xlabel ('Direction 1 (Largest Spread)');
ylabel ('Direction 2 (Next-Largest Spread)');
title (sprintf ('Principal Coordinates'));
arrow (scal*[1 0], "1");
arrow (scal*[0 1], "2");
for f=1:gmax  % loop over transcription factors f that INFLUENCE gene g
    arrow (scal*UgPrimeg([1 2], f), geneNames(f));
end
xlim ([-imageSize imageSize]);
ylim ([-imageSize imageSize]);


fprintf ("==================== READ THIS EXPLANATION! ======= \n");
fprintf ("In this example, the datapoints are distributed more or less as a flat ellipse.\n");
fprintf ("PCA identifies the two longest principal axis of this ellipse\n");
fprintf ("  (directions 1 and 2). \n");
fprintf ("\n");
return;

%======== COMPARE WITH STANDARD PCA?
% figure(3); clf; hold on;
% [coeff,score,latent] = pca (xkg);
% biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'H','K','G','N','B','C','T'});

 %============================================================
% NESTED UTILITY FUNCTIONS BELOW
%============================================================
   function arrow (dir, label)
      quiver (0, 0, dir(1), dir(2), 'LineWidth', 1.5);
      text (1.05*dir(1), 1.05*dir(2), label, 'FontSize', 16, 'FontWeight', 'Bold');      
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
end
