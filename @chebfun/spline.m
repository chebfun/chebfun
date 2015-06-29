function f = spline(x, y, d)
%SPLINE   CHEBFUN cubic spline data interpolation.
%   F = CHEBFUN.SPLINE(X, Y) returns a CHEBFUN F with domain [X(1), X(end)]
%   representing the cubic spline interpolant to the data values Y at the data
%   sites X. X must be a vector. If Y is a vector, then Y(j) is taken as the
%   value to be matched at X(j), hence Y must be of the same length as X  -- see
%   below for an exception to this. If Y is a matrix, then Y(:,j) is taken as
%   the value to be matched at X(j).
%
%   F = CHEBFUN.SPLINE(X, Y, D) is similar, but F is defined on the domain D.
%
%   Ordinarily, the not-a-knot end conditions are used. However, if Y contains
%   two more values than X has entries, then the first and last value in Y are
%   used as the end slopes for the cubic spline.
%
%   Example:
%   This generates a sine-like spline curve and samples it over a finer mesh:
%       x = 0:10;  y = sin(x);
%       f = chebfun.spline(x, y);
%       plot(x, y, 'o', f)
%
% See also INTERP1, PCHIP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Perhaps we shouldn't include the finals breaks in Not-A-Knot conditions?

if ( nargin < 3 )
    d = x([1,end]);
end

% Include breaks defined in the domain
breaks = unique([d(:) ; x(:)].');

% Number of intervals:
numInts = numel(breaks) - 1;

% Piecewise Chebyshev grid:
xx = chebpts(repmat(4, numInts, 1), breaks, 2);

% Forgive some transpose issues:
if ( ~any(length(x) + [0, 2] == size(y, 2)) )
    y = y.';
end

% Evaluate on the Chebyshev grid using built-in SPLINE:
yy = spline(x, y, xx);

% Orientate nicely:
if ( ~any(size(yy, 1) == 4*(length(x) + (-1:1))) )
    yy = yy.'; 
end

% Construct the CHEBFUN:
data = mat2cell(yy, repmat(4, numInts, 1), size(yy, 2));
f = chebfun(data, breaks, 'chebkind', 2);

% Restrict if needed:
if ( (d(1) > x(1)) || (d(end) < x(end)) )
    f = restrict(f, d);
end

end
