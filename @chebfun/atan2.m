function p = atan2(y, x, pref)
%ATAN2   Four quadrant inverse tangent of a CHEBFUN.
%   ATAN2(Y, X) is the four quadrant arctangent of the real parts of the CHEBFUN
%   objects X and Y.  -pi <= ATAN2(Y, X) <= pi.
%
%   ATAN2 is defined as:
%                  { atan(y/x),               x > 0
%                  { atan(y/x) + pi,  y >= 0, x < 0
%   atan2(y, x) =  { atan(y/x) - pi,  y < 0,  x < 0
%                  { pi/2,            y > 0,  x = 0
%                  { -pi/2,           y < 0,  x = 0
%                  { 0,               y = 0,  x = 0
%
% See also ATAN, ATAN2D.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Set up:
% Grab some preferences:
if ( nargin < 3)
    pref = chebfunpref();
end

if ( ~isreal(y) || ~isreal(x) )
    error('CHEBFUN:CHEBFUN:atan2:real', 'Inputs ust be real.');
end

numColsY = numColumns(y);
numColsX = numColumns(x);
if ( numColsY ~= numColsX )
    error('CHEBFUN:CHEBFUN:atan2:dim', 'Quasimatrix dimensions must agree.');
end

% Convert columns to a cell array:
y = num2cell(y);
x = num2cell(x);

% Loop over the columns:
for k = numel(y):-1:1
    p(k) = columnAtan(y{k}, x{k}, pref);
end

end

function p = columnAtan(y, x, pref)

% We'll need to extrapolate here:
pref.techPrefs.extrapolate = true;
tol = 2*eps*max(vscale(y), vscale(x));

% There's no reason why we shouldn't keep breaks in both x and y:
[x, y] = overlap(x, y);

%% Locate discontinuities:
% Find discontinuities resulting from y = 0:
ry = roots(y);

% But not where x > 0:
ry(feval(x, ry) > tol) = [];

% Deal with the case when y has zero FUNs:
rx = [];
for k = 1:numel(y.funs);
    % Here we look to see where x changes sign:
    if ( iszero(y.funs{k}) )
        rx = [ rx ; roots(x.funs{k}) ]; %#ok<AGROW>
    end
end

% Collect all the roots together:
r = unique([ ry ; rx ]);

% Remove near-by roots:
nr = numel(r) + 1;
while ( nr ~= numel(r) )
    nr = numel(r);
    r([diff(r) < tol ; false]) = [];
end
   
% Introduce new breaks at the computed roots if required:
if ( ~isempty(r) )
    newDom = union(x.domain, r.');
    x = restrict(x, newDom);
    y = restrict(y, newDom);
end

%% Compose:
% Do the composition:
p = compose(x, @(x, y) atan2(y, x), y, pref);

% Sort out the new point values:
if ( ~isempty(r) )
    % Set point values to zero if x(r) = y(r) = 0:
    idx = ((abs(x.pointValues) < tol) & (abs(y.pointValues) < tol));
    p.pointValues(idx) = 0;
end

end
