function f = pchip(x, y, d)
%PCHIP   CHEBFUN Cubic Hermite interpolating polynomial.
%   F = CHEBFUN.PCHIP(X, Y) returns a CHEBFUN F representing a certain
%   shape-preserving piecewise cubic Hermite interpolant to the values Y at the
%   sites X. X must be a vector. If Y is a vector, then Y(j) is taken as the
%   value to be matched at X(j), hence Y must be of the same length as X. If Y
%   is a matrix, then Y(:,j) is taken as the value to be matched at X(j).
%
%   F = CHEBFUN.PCHIP(X, Y, D) is similar, but F is defined on the domain D.
%
%  Example:
%    x = -3:3;
%    y = [-1 -1 -1 0 1 1 1];
%    plot(chebfun.pchip(x, y)); hold on, 
%    plot(chebfun.spline(x, y), '-r');
%    plot(x, y, 'or'), hold off
%    legend('pchip', 'spline')
%
% See also SPLINE, INTERP1.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    d = x([1,end]);
end

% Ensure x is a column vector:
if ( size(x, 1) == 1 )
    x = x.';
end

% Include breaks defined in the domain
breaks = unique([d(:) ; x(:)].');

% Number of intervals:
numInts = numel(breaks) - 1;

% Piecewise Chebyshev grid:
xx = chebpts(repmat(4, numInts, 1), breaks, 2);

% Forgive some transpose issues:
if ( ~any(size(y) == size(x))  )
    y = y.';
end

% Evaluate on the Chebyshev grid using built-in PCHIP:
if ( size(y, 2) == 1 ) % Note, PCHIP does weird things with orientation!
    yy = pchip(x, y, xx); 
else
    yy = pchip(x, y.', xx).'; 
end

% Construct the CHEBFUN:
data = mat2cell(yy, repmat(4, numInts, 1), size(yy, 2));
f = chebfun(data, breaks, 'chebkind', 2);

% Restrict if needed:
if ( (d(1) > x(1)) || (d(end) < x(end)) )
    f = restrict(f, d);
end

end
