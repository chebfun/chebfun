function h = heaviside(f)
%HEAVISIDE   Heaviside function of a CHEBFUN.
%   HEAVISIDE(F) returns a CHEBFUN which is -1 when F < 0, +1 when F > 0, and
%   0.5 when F == 0.
%
% See also DIRAC, SIGN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% HEAVISIDE() is basically a wrapper for SIGN().
h = .5*(sign(f) + 1);

end
