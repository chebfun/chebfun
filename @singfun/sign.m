function g = sign(f, pref)
%SIGN   Signum of a SINGFUN object.
%   SIGN(F) returns the sign of F, where F is a SINGFUN object with no 
%   roots in its domain. If F has roots, then SIGN(F) will return garbage
%   with no warning.
%
%   For the nonzero elements of complex F, SIGN(F) = F./ ABS(F).
%
% See also ABS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Extract the smoothPart of F:
g = f.smoothPart;
    
if ( isreal(f) )
    % Pick a random point:
    arbitraryPoint = 0.1273881594;

    g.values = sign(feval(g, arbitraryPoint));
    g.coeffs = g.values;
else
    % Call the SIGN function of the SMOOTHFUN:
    g = sign(g, pref);
end

end