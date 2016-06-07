function h = hscale(f)
%HSCALE   Horizontal scale of a CHEBFUN object.
%   HSCALE(F) returns the infinity norm of the domain of F if the domain of F is
%   bounded, and the value 1 if it is not.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Compute INF norm of the domain:
h = norm(f(1).domain([1, end]), inf);

% Unbounded domains are defined to have hscale = 1:
if ( isinf(h) )
    h = 1;
end

end
