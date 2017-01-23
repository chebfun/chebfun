function f = sign(f, pref)
%SIGN   Signum of a TRIGTECH object.
%   SIGN(F) returns the sign of F, where F is a TRIGTECH object with no 
%   roots in its domain. If F has roots, then SIGN(F) will return garbage
%   with no warning. 
%
%   For the nonzero elements of complex F, SIGN(F) = F ./ ABS(F).
%
% See also ABS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isreal(f) )
    % Evaluate at the two end points, and an arbitrary interior point:
    arbitraryPoint = 0.1273881594;
    fx = feval(f, [-1; arbitraryPoint; 1]);
    % Take the mean:
    meanfx = mean(fx, 1);
    % Compute the floor:
    f.coeffs = sign(meanfx);
    f.values = f.coeffs;
else
    f = compose(f, @(x) x./abs(x));
end

end
