function r = roots(f, varargin)
%ROOTS   Roots of a FUN in the interval [a,b].
%   ROOTS(F) returns the real roots of the FUN F in the interval [a,b].
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs (strings) and VALs may be any of the following:
%
%   ALL: 
%       [0] - Return only real-valued roots in [a,b].
%        1  - Return roots outside of [a,b] (including complex roots).
%
%   RECURSE:
%        0  - Compute roots without domain bisection (slower).
%       [1] - Bisect until length(F) < 50. (fast, but additional complex roots).
%
%   PRUNE:
%       [0]
%        1  - Prune 'spurious' complex roots if ALL == 1 and RECURSE == 0.
%
%   HSCALE:
%       [1] - Horizontal scale for adjusting relative tolerances.
%     double
%
%   If F is an array-valued BNDFUN then there is no reason to expect each column
%   to have the same number of roots. In order to return a useful output, the
%   roots of each column are computed and then padded with NaNs so that a matrix
%   may be returned. The columns of R = ROOTS(F) correspond to the columns of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    r = [];
    return
end

% Find the roots of the ONEFUN of f:
onefunRoots = roots(f.onefun, varargin{:});

% Map the roots found on [-1,1] to the interval [a,b]:
r = f.mapping.for(onefunRoots);

end