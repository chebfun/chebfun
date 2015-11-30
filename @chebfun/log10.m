function g = log10(f, pref)
%LOG10   Base 10 logarithm of a CHEBFUN.
%   LOG10(F) returns the base 10 logarithm of F. If F has an roots in its
%   domain, then the representation is likely to be inaccurate.
%
% See also LOG, LOG2, EXP, LOGM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

if ( nargin < 2 )
    pref = chebfunpref();
end

% Add breaks at the roots of f:
f = addBreaksAtRoots(f);

% Call COMPOSE():
g = compose(f, @(x) log10(x), [], pref);

end
