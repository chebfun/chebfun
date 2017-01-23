function [L, Lconst] = lebesgue(x, varargin)
%LEBESGUE   Lebesgue function for a set of interpolation points.
%   L = LEBESGUE(X), where X is a set of points in [-1, 1], returns the
%   Lebesgue function associated with polynomial interpolation in those points.
%
%   L = LEBESGUE(X, a, b) or LEBESGUE(X, [a, b]), where X is a set of points in
%   [a, b] returns the Lebesgue function associated with polynomial
%   interpolation in those points in that domain.
%
%   L = LEBESGUE(..., 'trig') does the same but for trigonometric interpolation
%   instead of polynomial interpolation.
%
%   [L, LCONST] = LEBESGUE(...) also returns the Lebesgue constant.
%
%   Example:
%     The following commands compare the Lebesgue functions and constants for 8
%     Chebyshev, Legendre, and equispaced points in [-1, 1]:
%
%     n = 8;
%     [L, c] = lebesgue(chebpts(n));
%     subplot(1, 3, 1),  plot(L),  title(['Chebyshev: ' num2str(c)])
%     grid on,  axis([-1 1 0 8])
%     [L, c] = lebesgue(legpts(n));
%     subplot(1, 3, 2),  plot(L),  title(['Legendre: ' num2str(c)])
%     grid on,  axis([-1 1 0 8])
%     [L, c] = lebesgue(linspace(-1, 1, n));
%     subplot(1, 3, 3),  plot(L),  title(['Equispaced: ' num2str(c)])
%     grid on,  axis([-1 1 0 8])

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Parse inputs.
[d, doTrig] = parseInputs(x, varargin{:});

% Construct the appropriate Lebesgue function.
if ( doTrig )
    L = trigLebesgue(x, d);
else
    L = polyLebesgue(x, d);
end

% Return the Lebesgue constant if asked.
if ( nargout == 2 )
    Lconst = norm(L, inf);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d, doTrig] = parseInputs(x, varargin)

% Parse inputs.  (NB:  We've already stripped out the X input.)
if ( nargin == 1 )     % LEBESGUE(X)
    d = [-1 1];
    doTrig = false;
elseif ( nargin == 2 ) % LEBESGUE(X, [A B]) or LEBESGUE(X, 'trig')
    if ( isnumeric(varargin{1}) )
        d = varargin{1};
        doTrig = false;
    elseif ( strcmpi(varargin{1}, 'trig') )
        d = [-1 1];
        doTrig = true;
    else
        error('CHEBFUN:lebesgue:parseInputs:badArg1', 'Invalid argument.');
    end
elseif ( nargin == 3 ) % LEBESGUE(X, A, B) or LEBESGUE(X, [A B], 'trig')
    if ( isnumeric(varargin{1}) && isnumeric(varargin{2}) )
        d = [varargin{1} varargin{2}];
        doTrig = false;
    elseif ( isnumeric(varargin{1}) && strcmpi(varargin{2}, 'trig') )
        d = varargin{1};
        doTrig = true;
    else
        error('CHEBFUN:lebesgue:parseInputs:badArg2', 'Invalid argument.');
    end
elseif ( nargin == 4 ) % LEBESGUE(X, A, B, 'trig')
    if ( isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
         strcmpi(varargin{3}, 'trig') )
        d = [varargin{1} varargin{2}];
        doTrig = true;
    else
        error('CHEBFUN:lebesgue:parseInputs:badArg3', 'Invalid argument.');
    end
else
    error('CHEBFUN:lebesgue:parseInputs:tooManyArgs', ...
          'Too many input arguments.');
end

if ( ~isequal(size(d), [1 2]) || (d(1) >= d(2)) )
    error('CHEBFUN:lebesgue:parseInputs:badDom', ...
          ['Domain input must be either two numbers or a row vector with ' ...
           'two elements, listed in ascending order.']);
end

if ( any(x < d(1) - 10*eps(d(1))) || any(x > d(2) + 10*eps(d(2))) )
    error('CHEBFUN:lebesgue:parseInputs:pointsOutsideDomain', ...
          sprintf(['Interpolation points must be inside domain ' ...
                   '([%.6f %.6f]).'], d(1), d(2)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = polyLebesgue(x, d)
%POLYLEBESGUE   Compute Lebesgue function for polynomial interpolation.
%  X - Interpolation nodes.
%  D - Interpolation domain.

% Compute the barycentric weights for the interpolation grid.
w = baryWeights(x);

% Set preferences.
pref = chebfunpref();
pref.techPrefs.sampleTest = false;
if ( isa(pref.tech(), 'chebtech') )
    % In between the interpolation nodes, the Lebesgue function is
    % guaranteed by definition to be a polynomial of degree at most
    % length(x) - 1.
    pref.techPrefs.fixedLength = length(x);
end

% Set breakpoints at the interpolation nodes.  (NB:  unique() returns the
% points in sorted order.)
dom = unique([x(:) ; d.']).';
L = chebfun(@(t) polyLebesgueFun(t, x(:), w), dom, pref);

if ( isa(pref.tech(), 'chebtech') )
    % Since we fixed the degree for chebtech-based representations, we must
    % call simplify manualy.
    L = simplify(L);
end

end

function L = polyLebesgueFun(t, x, w)
%POLYLEBESGUEFUN:  Evaluate Lebesgue function for poly. interp. grid at a point.
%  T - Evaluation points.
%  X - Interpolation nodes.
%  W - Barycentric weights.

% Based on barycentric formula.
L = ones(size(t)); % Note: L(x) = 1
mem = ismember(t, x);
for i = 1:numel(t)
    if ( ~mem(i) )
        xx = w./(t(i) - x);
        L(i) = sum(abs(xx))/abs(sum(xx));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = trigLebesgue(x, d)
%TRIGLEBESGUE   Compute Lebesgue function for trigonometric interpolation.
%  X - Interpolation nodes (odd-length vector).
%  D - Interpolation domain.

% NB:  This code only works if the number of points is odd.

% Map from the domain in d to [-pi, pi].
map = @(t) pi/(d(2) - d(1))*(2*t - d(1) - d(2));

% Map the points to [-pi, pi], since trigBaryWeights() assumes the points
% come from there.
xm  = map(x);

% Remove a periodic endpoint.  This will make it so that, e.g.,
% lebesgue(linspace(-pi, pi, 16), 'trig') does the "expected" thing instead
% of treating -pi and pi as two distinct points that are really, really
% close to each other.
if ( norm([xm(1), xm(end)] - [-pi, pi], Inf) < 2*pi*eps )
    xm(end) = [];
end

% We can't deal with even-length grids.
if ( mod(length(xm), 2) == 0 )
    error('CHEBFUN:lebesgue:trigLebesgue:evenLengthGrid', ...
          ['LEBESGUE for trigonometric interpolation requires an ' ...
           'odd-length grid.\n(If you supplied an odd-length grid, ' ...
           'perhaps it has a periodic endpoint?)']);
end

% Evaluate the barycentric weights.
w = trigBaryWeights(xm);

% Set breakpoints at the interpolation nodes.
dom = unique([x(:) ; d.']).';
L = chebfun(@(t) trigLebesgueFun(map(t), xm(:), w), dom);

end

function L = trigLebesgueFun(t, x, w)
%TRIGLEBESGUEFUN:  Evaluate Lebesgue func. for trig. interp. grid at a point.
%  T - Evaluation points.
%  X - Interpolation nodes (odd-length vector).
%  W - Trigonometric barycentric weights.

% NB:  This code only works if the number of points is odd.

% Based on the trigonometric barycentric formula.
L = ones(size(t)); % Note: L(x) = 1
mem = ismember(t, x);
for i = 1:numel(t)
    if ( ~mem(i) )
        xx = w./sin((t(i) - x)/2);
        L(i) = sum(abs(xx))/abs(sum(xx));
    end
end

end
