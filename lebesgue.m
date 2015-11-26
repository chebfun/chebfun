function [L, Lconst] = lebesgue(x, varargin)
%LEBESGUE   Lebesgue function for a set of interpolation points.
%   L = LEBESGUE(X), where X is a set of points in [-1, 1], returns the
%   Lebesgue function associated with polynomial interpolation in those points.
%
%   L = LEBESGUE(X, a, b) or LEBESGUE(X, [a, b]), where X is a set of points in
%   [a, b] returns the Lebesgue function associated with polynomial
%   interpolation in those points in that domain.
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

%  Copyright 2015 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Parse inputs.
if ( nargin == 1 )     % LEBESGUE(X)
    d = [-1 1];
elseif ( nargin == 2 ) % LEBESGUE(X, [A B])
    d = varargin{1};
elseif ( nargin == 3 ) % LEBESGUE(X, A, B)
    d = [varargin{1} varargin{2}];
else
    error('CHEBFUN:lebesgue:tooManyArgs', 'Too many input arguments.');
end

% Compute the barycentric weights for the interpolation grid.
w = baryWeights(x);

% Set preferences.
pref = chebfunpref();
pref.techPrefs.sampleTest = false;
if ( isa(pref.tech(), 'chebtech') )
    % In between the interpolation nodes, the Lebesgue function is guaranteed
    % by definition to be a polynomial of degree at most length(x) - 1.
    pref.techPrefs.fixedLength = length(x);
end

% Set breakpoints at the interpolation nodes.  (NB:  unique() returns the
% points in sorted order.)
dom = unique([x(:) ; d.']).';
L = chebfun(@(t) lebesgueFun(t, x(:), w), dom, pref);

if ( isa(pref.tech(), 'chebtech') )
    % Since we fixed the degree for chebtech-based representations, we must call
    % simplify manualy.
    L = simplify(L);
end

% Return the Lebesgue constant if asked.
if ( nargout == 2 )
    Lconst = norm(L, inf);
end

end

function L = lebesgueFun(t, x, w)
%LEBESGUEFUN:  Evaluate Lebesgue function for an interpolation grid at a point.
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
