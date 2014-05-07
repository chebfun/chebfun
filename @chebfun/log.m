function g = log(f, pref)
%LOG   Natural logarithm of a CHEBFUN.
%   LOG(F) returns the natural logarithm of F. If F has an roots in its domain,
%   then the representation is likely to be inaccurate.
%
% See also LOG1P, LOG2, LOG10, EXP, REALLOG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

if ( nargin < 2 )
    pref = chebfunpref();
end

% Add breaks at the roots of f:
f = addBreaksAtRoots(f);

% Call COMPOSE():
g = compose(f, @(x) log(x), [], pref);

end
