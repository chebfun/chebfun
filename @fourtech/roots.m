function out = roots(f, varargin)
%ROOTS   Roots of a FOURTECH in the interval [-pi,pi].
%   ROOTS(F) returns the real roots of the FOURTECH F in the interval [-pi,pi].
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs (strings) and VALs may be any of the following:
%
%   ALL: 
%       [0] - Return only real-valued roots in [-pi,pi].
%        1  - Return roots outside of [-pi,pi] (including complex roots).
%
%   RECURSE:
%        0  - Compute roots without interval subdivision (slower).
%       [1] - Subdivide until length(F) < 50. (causes additional complex roots).
%
%   PRUNE:
%       [0]
%        1  - Prune 'spurious' complex roots if ALL == 1 and RECURSE == 0.
%
%   If F is an array-valued FOURTECH then there is no reason to expect each
%   column to have the same number of roots. In order to return a useful output,
%   the roots of each column are computed and then padded with NaNs so that a
%   matrix may be returned. The columns of R = ROOTS(F) correspond to the
%   columns of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Simply create a chebfun of the FOURTECH and call min on it.
g = chebfun(@(x) f.feval(x),[-pi,pi]);
out = roots(g,varargin{:});

end
