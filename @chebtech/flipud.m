function f = flipud(f)
%FLIPUD   Flip/reverse a CHEBTECH object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Negate the odd coefficients:
f.coeffs(end-1:-2:1,:) = -f.coeffs(end-1:-2:1,:);
% Note: we use end-1:-2:1 rather than 1:2:end-1 as we only want to change the
% odd coefficients (and this avoids checking the length of the coeffs).

end
