function g = log2(f, pref)
%LOG2   Base 2 logarithm of a CHEBFUN.
%   LOG2(F) returns the base 2 logarithm of F. If F has an roots in its domain,
%   then the representation is likely to be inaccurate.
%
% See also LOG, LOG10, POW2, NEXTPOW2, REALMAX, REALMIN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

if ( nargin < 2 )
    pref = chebfunpref();
end

% Add breaks at the roots of f:
f = addBreaksAtRoots(f);

% Call COMPOSE():
g = compose(f, @(x) log2(x), [], pref);

end
