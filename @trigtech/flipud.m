function f = flipud(f)
%FLIPUD   Flip/reverse a TRIGTECH object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-PI,PI].

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Flip the values:
f.values = [ f.values(1,:); flipud(f.values(2:end,:)) ];

% Negate the odd coefficients:
f.coeffs = flipud(f.coeffs);
% Note: we use end-1:-2:1 rather than 1:2:end-1 as we only want to change the
% odd coefficients (and this avoids checking the length of the coeffs).

end
