function p = atan2(y, x, pref)
%ATAN2    Four quadrant inverse tangent of a chebfun.
%   ATAN2(Y, X) is the four quadrant arctangent of the real parts of the chebfun
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
%                { 0,               y == 0, x = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('CHEBFUN:atan2:notImplemented', ...
    'atan2() not yet implemented for chebfun objects.');

%% Grab some preferences:
%if ( nargin < 3)
%    pref = chebfun.pref();
%end
%
%% We'll need to extrapolate here:
%pref.misc.extrapolate = true;
%tol = 1e-6*max(x.vscale, y.vscale);
%
%% There's no reason why we shouldn't keep breaks in both x and y:
%[x, y] = overlap(x, y);
%
%% Find discontinuities resulting from y = 0:
%ry = roots(y);
%% But not where x = 0:
%ry(feval(x, ry ) > tol) = [];
%% Or where dy/dt = 0: (avoid introcducing a break at a double root)
%ry(abs(feval(diff(y), ry)) < tol) = [];
%
%% Deal with the case when y has zero funs:
%rx = [];
%for k = 1:numel(y.funs);
%    % Here we look to see where x changes sign:
%    if ( iszero(y.funs{k}) )
%        rx = [ rx ; x.funs{k} ]; %#ok<AGROW>
%    end
%end
%
%% Collect all the roots together:
%r = [ ry ; rx ];
%
%% Introduce new breaks at the computed roots if required:
%if ( ~isempty(r) )
%    newDom = union(x.domain, r.');
%    x = restrict(x, newDom);
%    y = restrict(y, newDom);
%end
%
%% Do the composition:
%p = compose(x, @(x, y) atan2(y, x), y, pref);
%
%% Sort out the new impulses:
%if ( ~isempty(r) )
%    % Find roots that are breakpoints and evaluate there:
%    [r, ignored, idx] = intersect(r', p.domain);
%    % Set these impulses to zero if x(r) = y(r) = 0:
%    z = abs(feval(x, r)) < tol & abs(feval(y, r)) < tol;
%    p.impulses(1, idx(z)) = 0;    % Zero where x = y = 0.
%    p.impulses(1, idx(~z)) = pi;  % Pi elsewhere.
%end

end
