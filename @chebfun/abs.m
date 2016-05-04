function F = abs(F, pref)
%ABS   Absolute value of a CHEBFUN.
%   ABS(F) is the absolute value of the CHEBFUN F.
%
% See also SIGN, ANGLE, UNWRAP, HYPOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case: (F is empty)
if ( isempty(F) )
    return
end

if ( nargin < 2 )
    pref = chebfunpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    F(k) = columnAbs(F(k), pref);
end

end

function g = columnAbs(f, pref)

% Add breaks at the appropriate roots of f:
g = addBreaksAtRoots(f, pref);

% Call ABS on each of the FUNs: (result will be smooth)
for k = 1:numel(g.funs)
    g.funs{k} = abs(g.funs{k});
end

% Take the absolute value of the point values at break points:
g.pointValues = abs(g.pointValues);

% [TODO]: Do we want to do this?
% [ignored, idx] = setdiff(f.domain, g.domain);
% g = merge(g, idx.', pref); 

% [TODO]: Do we want to do this?
g = simplify(g, pref);

% TODO: Should we always return a quasimatrix rather than overlap the roots?

end
