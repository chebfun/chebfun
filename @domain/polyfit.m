function f = polyfit(x, y, n, d)
%POLYFIT   Polyfit discrete data and return a CHEBFUN object.
%   F = POLYFIT(X, Y, N, D), where D is a DOMAIN object, returns a CHEBFUN F on
%   the domain D([1, end]) which corresponds to the polynomial of degree N that
%   fits the data (X, Y) in the least-squares sense. X should be a real-valued
%   column vector and Y should be a matrix with size(Y,1) = size(X,1).
%
%   Note DOMAIN/POLYFIT does not not support more than one output argument in
%   the way that MATLAB/POLYFIT does.
%
% See also CHEBFUN/POLYFIT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert domain to a double:
d = double(d);
d = d([1, end]);

% Align dimensions:
if ( size(x, 2) > 1 )
    if ( size(x, 1) > 1 )
        error('CHEBFUN:DOMAIN:polyfit:xIn', ...
            'X should be a real-valued column vector');
    end
    x = x.';
end
if ( size(y, 1) ~= size(x, 1) )
    if ( size(y, 2) ~= size(x, 1) )
        error('CHEBFUN:DOMAIN:polyfit:xIn', ...
            'X and Y vectors must be the same size.');
    end
    y = y.';
end

% Make Chebyshev-Vandermonde matrix:
T = chebpoly(0:n, d);
Tx = feval(T, x);
% Solve for coefficients (least squares)
c = Tx\y;
% Construct Chebfun:
f = chebfun(c, d, 'coeffs');

end