function f = flipud(f)
%FLIPUD   Flip/reverse a TRIGTECH object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Flip the values:
f.values = [ f.values(1,:); flipud(f.values(2:end,:)) ];

% Flip the coefficients taking into account where f is odd or even
if mod(size(f.coeffs,1),2)
    % Odd length is easy, just flip the coefficients
    f.coeffs = flipud(f.coeffs);
else
    % Even length requires keeping the first coefficient in place and 
    % flipping the remaining ones.  This follows since we interpret the 
    % first coefficient to correspond to the 1/2*cos(-N/2 x) mode.
    f.coeffs(1) = conj(f.coeffs(1));
    f.coeffs(2:end,:) = flipud(f.coeffs(2:end,:));
end

end
