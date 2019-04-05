function [] = plotPrincipalComponents (opts, xkg, yk, Tg, h)

fprintf ('\n\n\n\n===================================\n');

set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontWeight','Normal');  % Bold

%======== READ HARDWIRED FILES FOR NOW
numNuclei     = size (xntg,1);
numTimepoints = size (xntg,2);
numRegulators = size (xntg,3);
numDatas = numNuclei*numTimepoints;
xkg = reshape (xntg, [numDatas numRegulators]);
ykg = reshape (yntg, [numDatas  numRegulators]);
yk = ykg(:,g);

Tg = Thparams(1:numRegulators);
h = Thparams(end);

%======== CONSTRUCT 7-DIMENSIONAL ROTATION TO REVEAL STRUCTURE
% Store basis vectors in columns
basis = eye(numRegulators, numRegulators);
basis (:,1) = Tvec;
basis = gramSchmidtOrthogonalize (basis);
xorg = -h/norm(Tvec) * basis(:,1);
xkgRot = (xkg - repmat(xorg',[numDatas 1])) * basis;  % Transform all the points

% in R2018 can just write (xkg - xorg')

%======== CONSTRUCT 6-DIMENSIONAL ROTATION TO PRINCIPAL AXIS FRAME
% Direction 1 is fixed to be along T.
% From now on we only work with directions 2 to 7.
covMat = cov(xkgRot(2:end, 2:end));   % find 6x6 covariance matrix
[evecs,evals] = eig (covMat);
[evals,ind] = sort(diag(evals));   % ANNOYING - sort in descending order
evals = evals(ind);                % ANNOYING
evecs = evecs(:,ind);              % ANNOYING

secondRot = eye(numRegulators, numRegulators);
secondRot(2:end, 2:end) = evecs;

secondShift = -mean(xkgRot,1);     % find mean over all data points
secondShift(1) = 0;                % don't shift in direction 1

xkgRot = (xkgRot + repmat(secondShift,[numDatas 1])) * secondRot;



%======== BIPLOT SHOWING QUALITY OF CLASSIFICATION
hfig = figure(4);
set (hfig, 'pos',[0 600 300 300], 'Toolbar','None','MenuBar','None');
%print (sprintf('biplot%d',g), '-dpng', '-r300');

clf; hold on;
scatter (xkgRot(yk>0,1), xkgRot(yk>0,2), 20, [0 .6 0], 'Filled');
scatter (xkgRot(yk<0,1), xkgRot(yk<0,2), 20, [1 0 0]);
xlabel ('Principal direction 1 (switching boundary normal)');
ylabel ('Principal direction 2');
%title (sprintf ('Quality of Classification'));

scalfac = 100.;
for f=1:numRegulators
    quiver (0, 0, scalfac*basis(f,1), scalfac*basis(f,2), 'LineWidth', 2, 'Color', [0 0 0]);
    text (scalfac*basis(f,1), scalfac*basis(f,2), sprintf ('%s', opts.geneNames{f}), 'FontSize', 16, 'FontWeight', 'Bold');
end
xlim (1.2*[-scalfac scalfac]);
ylim (1.2*[-scalfac scalfac]);
plot ([0 0], [-1000 1000], '--', 'color', [.4 .4 .4]);



%======== DIAGNOSTICS
fprintf ('We have %d points in %d-dimensional gene expression space\n', numDatas, numRegulators);
fprintf ('We have on/off class data for each of these\n', numDatas, numRegulators);



% We need to collapse into 2D
% First let try Principal Component Analysis just for fun
% Let q1, q2, q3, etc, be the coordinates in the rotated basis

%figure(2); clf; hold on;
%[coeff,score,latent] = pca (xkg);
%biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'H','K','G','N','B','C','T'});






%=================================================
% Given a set of column vectors (c1,c2,...,cn),
% perform Gram-Schmidt orthogonalization:
%                                   c1 /= |c1|
% c2 = c2 - (c2.c1)c1;              c2 /= |c2|
% c3 = c3 - (c3.c1)c1 - (c3.c2)c2;  c3 /= |c3|
% and so on.
    function bNew = gramSchmidtOrthogonalize (b)
        assert (size(b,1) == size(b,2));
        n = size (b,1);
        for i=1:n        % vector to orthogonalize
            for j=1:i-1  % vector to orthogonalize against
                b(:,i) = b(:,i) - (b(:,i)' * b(:,j)) * b(:,j);
            end
            b(:,i) = b(:,i) ./ norm(b(:,i));  % normalize
        end
        bNew = b;
    end

end %function


function M=RandOrthMat(n, tol)
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

