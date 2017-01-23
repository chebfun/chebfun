function out = isSubset(A, B, tol)
%ISSUBSET   Test if A is a subset of B.
%   OUT = ISSUBSET(A, B, TOL) returns logical true if A is a subset of B up
%   to the tolerance TOL, where A and B are domains of a chebfun, chebfun2 
%   or chebfun3 (i.e., vectors of size 1x2, 1x4 or 1x6 with A(1) <= A(2), 
%   A(3) <= A(4), A(5) <= A(6)).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(A) )
    out = 1;
    return
end

if ( length(A) ~= length(B) )
    error('CHEBFUN:isSubset:size', ...
        'Domains must have the same number of entries.')
end

if ( ( A(1) < B(1) - tol ) || ( B(2) + tol < A(2) ) )
    out = 0;
    return
else
    % Recurse:
    out = isSubset(A(3:end), B(3:end), tol);
end

end