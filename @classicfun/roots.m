function r = roots(f, varargin)
%ROOTS   Roots of a CLASSICFUN in the interval [a,b].
%   ROOTS(F) returns the real roots of the CLASSICFUN F in the interval [a,b].
%
%   ROOTS(F, OPTIONS) modifies the default ROOTS properties, by passing the
%   OPTIONS to the rootfinding method of the ONEFUN of F.
%
%   If F is an array-valued CLASSICFUN then there is no reason to expect each column to
%   have the same number of roots. In order to return a useful output, the roots
%   of each column are computed and then padded with NaNs so that a matrix may
%   be returned. The columns of R = ROOTS(F) correspond to the columns of F.
%
% See also ONEFUN. [TODO]: Why may we also want to see ONEFUN?

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    r = [];
    return
end

% Find the roots of the ONEFUN of F:
onefunRoots = roots(f.onefun, varargin{:});

% Map the roots found on [-1,1] to the interval [a,b]:
r = f.mapping.For(onefunRoots);

end
