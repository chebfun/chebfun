function f = flipud(f)
%FLIPUD   Flip/reverse a CHEBTECH object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Negate the odd coefficients:
f.coeffs(2:2:end,:) = -f.coeffs(2:2:end,:);
% Note: the first odd coefficient is at index 2.

end
