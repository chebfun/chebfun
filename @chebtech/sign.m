function f = sign(f, pref)
%SIGN   Signum of a CHEBTECH object.
%   SIGN(F) returns the sign of F, where F is a CHEBTECH object with no 
%   roots in its domain.  If F has roots, then SIGN(F) will return garbage
%   with no warning. 
%
%   For the nonzero elements of complex F, SIGN(F) = F ./ ABS(F).
%
% See also ABS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( isreal(f) )
    arbitraryPoint = 0.1273881594;
    f.coeffs = sign(feval(f, arbitraryPoint));
%     f.coeffs = f.values;
    f.vscale = abs(f.coeffs);
else
    if ( nargin == 1 )
        pref = chebtech.techPref();
    end
    pref.extrapolate = 1;
    f = compose(f, @(x) x./abs(x), [], pref);
end

end
