function r = roots(f, varargin)
%ROOTS   Roots of a FUN in the interval [a,b].
%   ROOTS(F) returns the real roots of the FUN F in the interval [a,b].
%
%   ROOTS(F, OPTIONS) modifies the default ROOTS properties, by passing the
%   OPTIONS to the rootfinding method of the ONEFUN of F.
%
%   If F is an array-valued FUN then there is no reason to expect each column to
%   have the same number of roots. In order to return a useful output, the roots
%   of each column are computed and then padded with NaNs so that a matrix may
%   be returned. The columns of R = ROOTS(F) correspond to the columns of F.
%
%   See also ONEFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    r = [];
    return
end

% Find the roots of the ONEFUN of F:
onefunRoots = roots(f.onefun, varargin{:});

% Map the roots found on [-1,1] to the interval [a,b]:
r = f.mapping.for(onefunRoots);

%%
% Get rid of spurious roots which are caused by fast decay of function defined
% in an unbounded domain:

% TODO: NH: Move this to @UNBNDFUN/ROOTS()

% Set a threshold for the 'farfield':
farfield = 1e-1/eps;

ends = get(f, 'domain');

if ( isinf(ends(1)) )
    mask = ( r < -farfield );
    r(mask) = [];
end

if ( isinf(ends(2)) )
    mask = ( r > farfield );
    r(mask) = [];
end

end