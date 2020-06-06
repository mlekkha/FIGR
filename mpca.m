%============================================================
% Modified Principal Component Analysis for Binary-Classified Data
% Yen Lee Loh 2020-6-5
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
% USAGE:     [basis,xkg] = mpca (xkg, yk, Tg)
% ARGUMENTS: xkg(1:K,1:G)    coordinates of each data point
%            ykg(1:K)        class of each data point (+1 or -1)
%            Tg(1:G)         normal to classification hyperplane
% RETURNS:   basis(1:G,1:G)  principal-axis vectors
%            xkgRot(1:K,1:G) data coordinates in rotated & shifted frame
%
% The ROW vs COLumn stuff is ridiculously confusing!
%
%============================================================
function [finalRot, xkgRot] = mpca (xkg, yk, Tg, h)
fprintf ('\n============== mpca() =====================\n');

[kmax gmax] = size (xkg);

xgk = xkg';   % more natural if position vectors are column vectors
% after all, "basis" is going to be column vectors

%======== ROTATE TO ALIGN "X" AXIS WITH CLASSIFICATION HYPERPLANE NORMAL
% Direction 1 is fixed to be along Tg.
firstRot = eye (gmax, gmax);
firstRot (:,1) = Tg;
firstRot = gramSchmidtOrthogonalize (firstRot);
firstShf = -h/norm(Tg) * firstRot(:,1);
xgkRot = firstRot' * (xgk - firstShf);  
% Note the TRANSPOSE.
% The COORDINATES rotate in the OPPOSITE SENSE to the BASIS vectors.

fprintf ("Intermediate basis (1st column is Tg):\n"); disp (firstRot);

%======== CONSTRUCT 6-DIMENSIONAL ROTATION TO PRINCIPAL AXIS FRAME
% Rotate within the subspace spanned by directions 2,3,4,...,G.
covMat = cov (xgkRot'); % find covariance matrix
fprintf ("Covariance matrix:\n"); disp (covMat);

covMat = covMat(2:end,2:end);
fprintf ("Covariance matrix for directions 2...G:\n"); disp (covMat);

[evecs,evals] = eig (covMat);  % find eigenvalues and eigenvectors
evals = diag (evals);                        % of (G-1)x(G-1) lower right part of cov

[~,ind] = sort (evals);  % find indices that sort evals in ascending order
ind = flip (ind);
evals = evals(ind);      % sort evals in descending order of evals
evecs = evecs(:,ind);    % sort evecs in descending order of evals

fprintf ("Eigenvals in decreasing order:\n"); disp (evals'); 
fprintf ("Eigenvecs (columns) in decreasing order of evals:\n"); disp (evecs);

secondRot = eye (gmax, gmax);
secondRot(2:end,2:end) = evecs;
secondShf = mean(xgkRot,2);      % find mean over all data points
secondShf(1) = 0;                % don't shift in direction 1

fprintf ("Rotation within subspace spanned by dirs 2...G:\n"); disp (secondRot);

xgkRot = secondRot' * (xgkRot - secondShf);
finalRot = firstRot * secondRot;       % basically (secondRot'*firstRot')'

fprintf ("Row vectors below are modified principal axes:\n"); disp (finalRot);

xkgRot = xgkRot';  % convert back to correct form for returnvalue
fprintf ("============== end mpca() =====================\n");
return  % return [finalRot, xkgRot]

% If we wished, we could have concatenated the two transformations into
% a single transformation of the form
% xgkRot = finalRot' * (xgkRot - finalOrg);
%
%
%



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
