function C = tolUnique(A, tol)
%TOLUNIQUE   UNIQUE with a tolerance for checking floating-point equality.
%   C = TOLUNIQUE(A, TOL) for a real vector A is the same as UNIQUE(A) except
%   that two values are deemed "equal" if they are within an absolute
%   difference TOL of one another.  If this happens, the mean of the two values
%   is placed in C.
%
%   C = TOLUNIQUE(A) uses a default tolerance of 100*EPS*MAX(NORM(A, Inf)).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    tol = 100*eps*max(norm(A, Inf));
end

% TODO:  This might not do the right thing if there are multiple values that
% are all close to one another, say, 2, 2 + 2*eps, 2 + 4*eps, and 2 + 6*eps,
% and tol = 3*eps.  This code will merge all of those into one value:  2 + eps
% (the average of 2 and 2 + 2*eps), and one can argue that behavior is wrong.
% We can resolve this if this ever becomes an issue in practice.
C = unique(A);
ind = find(diff(C) < tol);
C(ind) = (C(ind) + C(ind+1))/2;
C(ind+1) = [];

end
