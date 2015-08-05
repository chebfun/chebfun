function g = reallog(f, pref)
%REALLOG   Real logarthm of a CHEBFUN.
%   REALLOG(F) is the natural logarithm of the elements of F. An error is
%   produced if F is negative. If F has an roots in its domain, then the
%   representation is likely to be inaccurate.
%
% See also LOG, LOG2, LOG10, EXP, REALPOW, REALSQRT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

% Check for complex CHEBFUN objects:
if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:reallog:complex', ...
        'Reallog produced complex result.');
end

% Obtain preferences:
if ( nargin < 2 )
    pref = chebfunpref();
end

% Add breaks at the roots of f:
f = addBreaksAtRoots(f);

% Call COMPOSE():
g = compose(f, @(x) reallog(x), [], pref);

end
