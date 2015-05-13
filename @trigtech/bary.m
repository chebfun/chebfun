function fx = bary(x, fvals, xk, wk)
%BARY   Trigonometric barycentric interpolation formula.
%   TRIGTECH.BARY(X, FVALS, XK, VK) uses the 2nd form trigonometric barycentric
%   formula with weights VK to evaluate an interpolant of the data {XK,
%   FVALS(:,k)} at the points X. Note that XK must be equally-spaced.
%   Furthermore, XK and VK should be column vectors, and FVALS, XK, and VK
%   should have the same length.
%
%   BARY(X, FVALS) assumes XK are equispaced in [-1, 1) (i.e., as return by
%   TRIGTECH.TRIGPTS()) and that VK are the corresponding barycentric weights
%   for trigonomentric-polynomial interpolation (i.e, as returned by
%   TRIGTECH.BARYWTS().).

% See also TRIG.TRIGPTS and TRIG.BARYWTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Length of interpolant:
n = length(fvals);

% Store dimension of input and reshape to column vector:
sizex = size(x);
x = x(:);

% Default to equigrid:
if ( nargin < 3 )
    xk = trigtech.trigpts(n);
end
% Default weights to trig-poly on equigrid:
if ( nargin < 4 )
    wk = trigtech.barywts(n);
end

% Scale to an interval of length 2*pi:
dom = xk(end) + xk(2) - 2*xk(1);
scl = 2*pi/dom;

% Deal with odd / even points (and include sclaing):
if ( mod(n, 2) )
    F = @(x) sin(scl*x); % Odd number of points
else
    F = @(x) tan(scl*x); % Even number of points 
end   

% Intialise:
numer = 0; denom = 0;
% Loop over the interpolation nodes:
for k = 1:n
    temp = wk(k) ./ F((x-xk(k))/2);
    numer = numer + temp*fvals(k);
    denom = denom + temp;
end
fx = numer./denom;

% Try to clean up NaNs:
for k = find(isnan(fx(:,1)))'       % (Transpose as Matlab loops over columns)
    index = find(x(k) == xk, 1);    % Find the corresponding node
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
    end
end

% Reshape output:
fx = reshape(fx, sizex);

end