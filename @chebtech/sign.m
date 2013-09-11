function f = sign(f, pref)
%SIGN   Signum of a CHEBTECH object.
%   SIGN(F) returns the absolute value of F, where F is a CHEBTECH object with
%   no roots in F.domain. If ~isempty(roots(F)) then SIGN(F) will return garbage
%   with no warning.
%
%   For the nonzero elements of complex F, sign(F) = F ./ ABS(F).
%
% See also ABS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isreal(f) )
    arbitraryPoint = 0.1273881594;
    f.values = sign(feval(f, arbitraryPoint));
    f.coeffs = f.values;
    f.vscale = abs(f.values);
else
    if ( nargin == 1 )
        pref = chebtech.pref;
    end
    pref.chebtech.extrapolate = 1;
    f = compose(f, @(x) x./abs(x), [], pref);
end

end