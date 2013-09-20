function p = unwrap(p, jumptol)
%UNWRAP Unwrap CHEBFUN phase angle.
%   UNWRAP(P) unwraps radian phases P by changing absolute jumps greater than or
%   equal to pi to their 2*pi complement. It unwraps along the continuous
%   dimension of P and leaves the first FUN along this dimension unchanged.
%
%   UNWRAP(P, TOL) uses a jump tolerance of TOL rather than the default TOL =
%   pi.
%
%   See also ABS, ANGLE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trivial case
if ( isempty(p) || numel(p.funs) == 1 )
    return
end

if ( size(f.funs{1}, 2) > 1 )
    error('CHEBFUN:unwrap:array', ...
        'UNWRAP() does not support array-valued CHEBFUN objects.');
end

% Get the shifts. By default these are pi.
if ( nargin == 1 )
    jumptol = pi;
end

% Choose a tolerance:
vs = get(p, 'vscale');
el = get(p, 'epslevel');
tol = 100*el*vs;

% Store data about the impulses for later:
lvals = feval(p, p.domain, 'left');
rvals = feval(p, p.domain, 'right');
idxl = abs(p.impulses(:,:,1) - lvals) < 100*tol;
idxr = abs(p.impulses(:,:,1) - rvals) < 100*tol;

% Find the jumps:
idx1 = mymod(abs(lvals - rvals), 2*jumptol) <  tol;
idx2 = abs(lvals - rvals) > tol;
idx = idx1 & idx2;
% Scale and shift:
scale = round((lvals - rvals) / (2*jumptol));
shift = cumsum(idx.*scale*2*jumptol);
% Update the FUNs:
for j = 2:numel(p.funs)
    p.funs{j} = p.funs{j} + shift(j);
end

% Update the impulses:
p.impulses(1,idxl) = feval(p, p.domain(idxl), 'left');
p.impulses(1,idxr) = feval(p, p.domain(idxr), 'right');

% Merge to tidy up unneeded breakpoints:
p = merge(p, find(idx));

function m = mymod(f, g)
    m = min([abs(mod(f,g)) ; abs(mod(f,-g)) ; abs(mod(-f,g)) ; abs(mod(-f,-g))]);
end

end