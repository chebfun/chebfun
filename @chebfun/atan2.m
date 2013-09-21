function p = atan2(y, x, pref)
%ATAN2    Four quadrant inverse tangent of a CHEBFUN.
%   ATAN2(Y, X) is the four quadrant arctangent of the real parts of the CHEBFUN
%   objects X and Y.  -pi <= ATAN2(Y,X) <= pi.
%
% See also ATAN, ATAN2D.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATAN2 is defined as:
%                { atan(y/x),               x > 0
%                { atan(y/x) + pi,  y >= 0, x < 0
% atan2(y, x) =  { atan(y/x) - pi,  y < 0,  x < 0
%                { pi/2,            y > 0,  x = 0
%                { -pi/2,           y < 0,  x = 0
%                { 0,               y = 0,  x = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error('CHEBFUN:atan2:notImplemented', ...
%     'atan2() not yet implemented for CHEBFUN objects.');

% Grab some preferences:
if ( nargin < 3)
    pref = chebfun.pref();
end

% We'll need to extrapolate here:
pref.chebfun.extrapolate = true;
tol = max(epslevel(y)*vscale(y), epslevel(x)*vscale(x));

% There's no reason why we shouldn't keep breaks in both x and y:
[x, y] = overlap(x, y);

% Get the number of columns:
numCols = min(size(x));

%%
% Find discontinuities resulting from y = 0:
ry = roots(y);

% But not where x = 0:
xr = feval(x, ry);
% Deal with array-valued CHEBFUN objects:
xr = xr(:, 1:numCols);
% Remove unwanted roots:
ry(xr > tol) = [];

% Or where dy/dt = 0: (avoid introcducing a break at a double root)
ypr = feval(diff(y), ry);
% Deal with array-valued CHEBFUN objects:
ypr = ypr(:, 1:numCols);
% Remove unwanted roots:
ry(abs(ypr) < tol) = [];

% Deal with array-valued CHEBFUN objects:
ry = unique(ry(:));
ry(isnan(ry)) = [];

%%
% Deal with the case when y has zero FUNs:
rx = [];
for k = 1:numel(y.funs);
    % Here we look to see where x changes sign:
    if ( iszero(y.funs{k}) )
        rx = [ rx ; roots(x.funs{k}) ]; %#ok<AGROW>
    end
end
% Deal with array-valued CHEBFUN objects:
rx = unique(rx(:));
rx(isnan(rx)) = [];

%%
% Collect all the roots together:
r = [ ry ; rx ];

% Introduce new breaks at the computed roots if required:
if ( ~isempty(r) )
    newDom = union(x.domain, r.');
    x = restrict(x, newDom);
    y = restrict(y, newDom);
end

%%
% Do the composition:
p = compose(x, @(x, y) atan2(y, x), y, pref);

% Sort out the new impulses:
if ( ~isempty(r) )
    % Find roots that are breakpoints and evaluate there:
    [r, ignored, idx] = intersect(r', p.domain);
    idx = repmat(idx, 1, numCols);
    % Set these impulses to zero if x(r) = y(r) = 0:
    z = abs(feval(x, r)) < tol & abs(feval(y, r)) < tol;
    
    idx1 = idx(z);
    if ( ~isempty(idx1) )
        % Zero where x = y = 0.
        for k = 1:numCols
            p.impulses(idx1(:,k),k,1) = 0;
        end
    end
    
    idx2 = idx(~z);
    if ( ~isempty(idx2) )
        % Pi elsewhere.
        for k = 1:numCols
            p.impulses(idx2(:,k),k,1) = pi;
        end
    end
end

end
