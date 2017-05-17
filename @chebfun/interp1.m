function p = interp1(x, y, method, dom)
%INTERP1   CHEBFUN polynomial interpolant at any distribution of points.
%   P = CHEBFUN.INTERP1(X, Y), where X and Y are vectors, returns the CHEBFUN P
%   defined on the domain [X(1), X(end)] corresponding to the polynomial
%   interpolant through the data Y(j) at points X(j).
%
%   If Y is a matrix with more than one column then Y(:,j) is taken as the
%   value to be matched at X(j) and P is an array-valued CHEBFUN with each
%   column corresponding to the appropriate interpolant.
%
%   EXAMPLE: The following commands plot the interpolant in 11 equispaced points
%   on [-1, 1] through the famous Runge function:
%       d = [-1, 1];
%       ff = @(x) 1./(1+25*x.^2);
%       x = linspace(d(1), d(2), 11);
%       p = chebfun.interp1(x, ff(x))
%       f = chebfun(ff, d);
%       plot(f, 'k', p, 'r-'), hold on, plot(x, ff(x), 'r.'), grid on
%
%   P = CHEBFUN.INTERP1(X, Y, METHOD) specifies alternatative interpolation
%   methods.  The default is as described above. (Use an empty matrix []
%   to specify the default.) Available methods are:
%       'linear'   - linear interpolation
%       'spline'   - piecewise cubic spline interpolation (SPLINE)
%       'pchip'    - shape-preserving piecewise cubic interpolation
%       'cubic'    - same as 'pchip'
%       'poly'     - polynomial interpolation, as described above
%       'trig'     - trigonometric polynomial interpolation, as above
%       'periodic' - same as the 'trig' option.
%   
%   For the trigonometric case, if the endpoints of the domain coincide
%   with the first and the last interpolation points, the average of the 
%   corresponding function values is interpolated.
% 
%   P = CHEBFUN.INTERP1(X, Y, METHOD, DOM) restricts the result P to the domain
%   DOM.
%
% See also SPLINE, PCHIP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin == 2 )
    method = [];
    dom = [];
end

if ( nargin == 3 )
    if ( ischar(method) )
        dom = [];
    else
        dom = method;
        method = [];
    end
end
dom = double(dom);

if ( isempty(method) )
    % Default to poly:
    method = 'poly';
end

% Ensure x and y are both column vectors:
if ( size(x, 1) == 1 )
    x = x.';
end

if ( size(y, 1) == 1 )
    y = y.';
end

% Sort x and also y if y is not a chebfun:
[x, idx] = sort(x);
if ( isa(y, 'chebfun') )
    y = feval(y, x);
else
    y = y(idx,:);
end

% Set default domain if none was supplied.
if ( isempty(dom) )
    dom = [x(1) x(end)];
end

switch method
    case 'poly'
        p = interp1Poly(x, y, dom);
    case {'trig', 'periodic'}
        p = interp1Trig(x, y, dom);
    case 'spline'
        p = chebfun.spline(x, y, dom);
    case {'pchip', 'cubic'}
        p = chebfun.pchip(x, y, dom);
    case 'linear'
        p = interp1Linear(x, y, dom);
    otherwise
        error('CHEBFUN:CHEBFUN:interp1:method', 'Unknown method ''%s''', ...
            method);
end

end

function p = interp1Poly(x, y, dom)
% Polynomial interpolation

% Compute barycentric weights for these points:
w = baryWeights(x);
% Define the interpolant using CHEBTECH.BARY():
f = @(z) bary(z, y, x, w);
% Construct a CHEBFUN:
p = chebfun(f, dom, length(x));

end

function p = interp1Trig(x, y, dom)
% Trigonometric interpolation

% Remove a periodic end-point and interpolate the average:
if ( norm([x(1), x(end)] - dom)/diff(dom) < 100*eps )
    x(end) = [];
    y(1, :) = (y(1, :) + y(end, :))/2;
    y(end, :) = [];
end

n = length(x);
% Evaluate the interpolant on n equally spaced points using 
% the trigonometric barycentric formula:
xx = trigpts(n, dom);
fx = trigBary(xx, y, x, dom);
% Construct a CHEBFUN:
p = chebfun(fx, dom, 'trig');
end


function p = interp1Linear(x, y, dom)
% Linear interpolation

% Include breaks defined in the domain
breaks = unique([dom(:) ; x(:)].');

% Number of intervals:
numInts = numel(breaks) - 1;

% Piecewise Chebyshev grid:
xx = chebpts(repmat(2, numInts, 1), breaks, 2).';

% Evaluate on the Chebyshev grid using built-in INTERP1:
yy = interp1(x, y, xx.', 'linear');

% Construct the CHEBFUN:
data = mat2cell(yy, repmat(2, numInts, 1), size(yy, 2));
p = chebfun(data, breaks, 'chebkind', 2);

% Restrict if needed:
if ( (dom(1) > x(1)) || (dom(end) < x(end)) )
    p = restrict(p, dom);
end

end
