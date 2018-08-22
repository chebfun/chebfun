function b = iszero( f )
%ISZERO   Check if a BALLFUN is identically zero on its domain.
%   OUT = ISZERO( F ) return 1 if the BALLFUN is exactly the zero function, and
%   0 otherwise.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if f = 0 at machine precision eps
eps = 1e-10;
% Test on the coeffs
b = max(max(max(abs(f.coeffs))))<eps;
% Test on the values
%b = max(max(max(abs(ballfun.coeffs3vals(f.coeffs)))))<eps;
end
