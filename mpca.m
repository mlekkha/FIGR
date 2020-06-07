%============================================================
% Modified Principal Component Analysis for Binary-Classified Data
% Yen Lee Loh 2020-6-6
%
% Given a set of datapoints in G-dimensional space
% classified into two groups,
% this routine constructs a set of G orthogonal unit vectors.
% The first vector (the "major principal direction" or "major axis") is
%   the normal to the boundary (hyperplane) that best separates the two classes.
% The subsequent (G-1) vectors are the eigenvectors of the covariance matrix
%   of the projections transverse to the major axis.
% In summary, the MPCA algorithm resolves datapoints
%   primarily in the "direction of classification",
%   and secondarily in the "directions of most spread".
%
% USAGE:     [UggP,xkgP] = mpca (xkg, yk, Tg, h)
% ARGUMENTS: xkg(1:K,1:G)    coordinates x_{kg} of each data point
%            yk(1:K)         class of each data point, y_k = +1 or -1
%            Tg(1:G)         normal T_g to classification hyperplane
%            h               offset parameter
% RETURNS:   UggP(1:G,1:G)   principal-axis vectors U_{gg'}
%            xkgP(1:K,1:G)   coordinates x_{kg'} in principal frame
%
% The algorithm works by constructing a coordinate transformation from
% x_g to x_{g''}, where the 1'' axis is aligned with the decision boundary normal,
% and then a second coordinate transformation to x_{g'}, which
% is a PCA in the reduced subspace spanned by the 2''...G'' axes.
%============================================================
function [UggP, xkgP] = mpca (xkg, yk, Tg, h)
fprintf ('\n============== mpca() =====================\n');

[kmax gmax] = size (xkg);

xgk = xkg';   % more natural if position vectors are column vectors
% after all, "basis" is going to be column vectors

%======== ROTATE TO ALIGN "X" AXIS WITH CLASSIFICATION HYPERPLANE NORMAL
% Direction 1 is fixed to be along Tg.
% Go from the unprimed frame to the double-primed (PP) frame.
UggPP = eye (gmax, gmax);
UggPP (:,1) = Tg;
UggPP = gramSchmidtOrthogonalize (UggPP);  % U_{g,g''}
UgPPg = inv (UggPP);
xgkMean = -h/norm(Tg) * UggPP(:,1);
xgPPk = UgPPg * (xgk - xgkMean);

fprintf ("Intermediate basis (1st column is Tg):\n"); disp (UggPP);

%======== CONSTRUCT 6-DIMENSIONAL ROTATION TO PRINCIPAL AXIS FRAME
% Rotate within the subspace spanned by directions 2,3,4,...,G.
Cgg = cov (xgPPk'); % find covariance matrix
fprintf ("Covariance matrix:\n"); disp (Cgg);

Cgg = Cgg(2:end,2:end);
fprintf ("Covariance matrix for directions 2...G:\n"); disp (Cgg);

[evecs,evals] = eig (Cgg);  % find eigenvalues and eigenvectors
evals = diag (evals);       % of (G-1)x(G-1) lower right part of C_{g_1,g_2}
[~,ind] = sort (evals);  % find indices that sort evals in ascending order
ind = flip (ind);
evals = evals(ind);      % sort evals in descending order of evals
evecs = evecs(:,ind);    % sort evecs in descending order of evals
evecs(:,2) = sign(det(evecs))*evecs(:,2); % hack to make det U=+1

fprintf ("Eigenvals in decreasing order:\n"); disp (evals'); 
fprintf ("Eigenvecs (columns) in decreasing order of evals:\n"); disp (evecs);

UgPPgP = eye (gmax,gmax);
UgPPgP(2:end,2:end) = evecs;
xgPPkMean = mean(xgPPk,2);      % find mean over all data points
xgPPkMean(1) = 0;                % don't shift in direction 1

fprintf ("Rotation within subspace spanned by dirs 2...G:\n"); disp (UgPPgP);

xgPk = UgPPgP' * (xgPPk - xgPPkMean);
UggP = UggPP * UgPPgP;       % U_{gg'} = U_{gg''} U_{g''g'}

fprintf ("Row vectors below are modified principal axes:\n"); disp (UggP);

xkgP = xgPk';  % convert back to correct form for returnvalue
fprintf ("============== end mpca() =====================\n");
return  % return [UggP, xkgP]

% If we wished, we could have concatenated the two transformations into
% a single transformation of the form
% xgkRot = finalRot' * (xgkRot - finalOrg);


    function bNew = gramSchmidtOrthogonalize (b)
        %=================================================
        % Given a set of column vectors (c1,c2,...,cn),
        % perform Gram-Schmidt orthogonalization as follows:
        %                                   c1 /= |c1|
        % c2 = c2 - (c2.c1)c1;              c2 /= |c2|
        % c3 = c3 - (c3.c1)c1 - (c3.c2)c2;  c3 /= |c3|
        % 
        % and so on, to obtain a set of n orthonormal vectors.        
        %=================================================
        assert (size(b,1) == size(b,2));
        n = size (b,1);
        for i=1:n        % vector to orthogonalize
            for j=1:i-1  % vector to orthogonalize against
                b(:,i) = b(:,i) - (b(:,i)' * b(:,j)) * b(:,j);
            end
            b(:,i) = b(:,i) ./ norm(b(:,i));  % normalize
        end
        bNew = b;     
        %=================================================
        % FOR DEBUGGING PURPOSES: Verify orthogonality
        %=================================================
        assert ( norm(bNew*bNew' - eye(n)) < 1e-14 );
    end
end
