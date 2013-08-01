function h = hscale(f)
%HSCALE    Horizontal scale of a CHEBFUN object.
%   HSCALE(F) returns the infinity norm of the domain of F if the domain of F is
%   bounded, and the value 1 if it is not. The horizontal scales of the
%   piecewise compnents of F are returned by get(F, 'hscale-local');

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Compute INF norm of the domain:s
h = norm(f.domain, inf);

% Unbounded domains are defined to have hscale = 1:
if ( isinf(h) )
    h = 1;
end

end