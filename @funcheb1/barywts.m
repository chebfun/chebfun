function v = barywts(n)
%BARYWTS % barycentric weights for Chebyshev points of 1st kind.
%   BARYWTS(N) returns the N barycentric weights for polynomial interpolation on
%   a Chebyshev grid of the 1st kind.
%
% See also FUNCHEB1.BARY.m, FUNCHEB1.CHEBPTS. 

% REFERENCE
%
% [1] Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange
% Interpolation", SIAM Review 46 (2004), pp 501-517.

v = sin( (2*(n-1:-1:0)+1)*pi/(2*n) ).';
v(2:2:end) = -v(2:2:end);
