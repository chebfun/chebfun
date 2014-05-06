function g = log1p(f, pref)
%LOG1P   Compute LOG(1+F) of a CHEBFUN F accurately.
%   LOG1P(F) computes LOG(1+F) without computing 1+F for small F. If F+1 has an
%   roots in its domain, then the representation is likely to be inaccurate.
%
% See also LOG, EXPM1.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

if ( nargin < 2 )
    pref = chebfunpref();
end

% Add breaks at the roots of f + 1:
f = addBreaks(f, getRootsForBreaks(f + 1));

% Call COMPOSE():
g = compose(f, @(x) log1p(x), [], pref);

end
