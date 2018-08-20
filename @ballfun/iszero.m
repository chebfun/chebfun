function b = iszero(f)
% ISZERO Test the nullity of a BALLFUN function up to machine
% precision
%   ISZERO(f) is the boolean |f| < 1e-10

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if f = 0 at machine precision eps
eps = 1e-10;
% Test on the coeffs
b = max(max(max(abs(f.coeffs))))<eps;
% Test on the values
%b = max(max(max(abs(ballfun.coeffs3vals(f.coeffs)))))<eps;
end
