function [x, idxV, idxW] = numIntersect(V, W, tol)
%NUMINTERSECT   Find the intersection of V and W with a relative tolerance.
%   The code does not always return unique elements.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 || isempty(tol) )
    pref = chebfunpref();
    tol = pref.deltaPrefs.proximityTol;
end
    
% Make sure V and W are vectors
V = V(:);
W = W(:);

% Intersection can not have more elements than either V or W:
m = min(length(V), length(W));
x = zeros( 1, m);
idxV = zeros(1, m);
idxW = zeros(1, m);

k = 0;
for i = 1:length(V)
    absVi = abs(V(i));
    for j = 1:length(W)
        absWj = abs(W(j));
        if ( V(i) == W(j) ) % Also handles when both are zero.
            % If they are exactly equal:
            x(k + 1) = V(i);
            idxV(k + 1) = i;
            idxW(k + 1) = j;
            k = k + 1;
        elseif ( abs(V(i) - W(j))/max(absWj, absVi) < tol )
            % If they are numerically equal:
            x(k + 1) = V(i);
            idxV(k + 1) = i;
            idxW(k + 1) = j;
            k = k + 1;
        end
    end
end

% Return valid entries only: 
if ( k == 0 )
    x = [];
    idxV = [];
    idxW = [];
else
    x = x(1:k);
    idxV = idxV(1:k);
    idxW = idxW(1:k);
end

end
