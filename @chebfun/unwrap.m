function p = unwrap(p, jumptol, pref)
%UNWRAP Unwrap chebfun phase angle.
%   UNWRAP(P) unwraps radian phases P by changing absolute jumps greater than or
%   equal to pi to their 2*pi complement. It unwraps along the continuous
%   dimension of P and leaves the first fun along this dimension unchanged.
%
%   UNWRAP(P, TOL) uses a jump tolerance of TOL rather than the default TOL =
%   pi.
%
%   See also UNWRAP, CHEBFUN/ABS, CHEBFUN/ANGLE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Trivial case
if ( numel(p.funs) == 1 )
    return
end

% Get the shifts. By default these are pi.
if ( nargin == 1 )
    jumptol = pi;
    pref = chebfun.pref();
elseif ( nargin == 1 )
    if isstruct(jumptol)
        pref = jumptol;
        jumptol = pi;
    else
        pref = chebfun.pref();
    end
end

% Store data about the imps for later:
vs = get(p, 'vscale');
vs = max([vs{:}]);
tol = 100*pref.chebfun.eps*vs;
lvals = feval(p, p.domain, 'left');
rvals = feval(p, p.domain, 'right');
idxl = abs(p.impulses(1,:) - lvals) < 100*tol;
idxr = abs(p.impulses(1,:) - rvals) < 100*tol;

idx1 = mymod(abs(lvals - rvals), 2*jumptol) <  tol;
idx2 = abs(lvals - rvals) > tol;
idx = idx1 & idx2;
scale = round((lvals - rvals) / (2*jumptol));
shift = cumsum(idx.*scale*2*jumptol);
for j = 2:numel(p.funs)
    p.funs{j} = p.funs{j} + shift(j);
end

% Update the imps
p.impulses(1,idxl) = feval(p, p.domain(idxl), 'left');
p.impulses(1,idxr) = feval(p, p.domain(idxr), 'right');

% Merge to tidy up unneeded breakpoints;
p = merge(p, find(idx));

function m = mymod(f,g)
m = min( [abs(mod(f,g)) ; abs(mod(f,-g)) ; abs(mod(-f,g)) ; abs(mod(-f,-g)) ]);
end

end