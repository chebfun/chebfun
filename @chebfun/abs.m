function g = abs(f)
%ABS   Absolute value of a CHEBFUN.
%   ABS(F) is the absolute value of the CHEBFUN F.
%
% See also SIGN, ANGLE, UNWRAP, HYPOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trivial case: (f is empty)
if ( isempty(f) )
    return
end

% Add breaks at the appropriate roots of f:
g = addBreaksAtRoots(f);

% Call ABS on each of the FUNs: (result will be smooth)
for k = 1:numel(g.funs)
    g.funs{k} = abs(g.funs{k});
end

% Take the absolute value of the point values at break points:
g.pointValues = abs(g.pointValues);

% [TODO]: Do we want to do this?
% [ignored, idx] = setdiff(f.domain, g.domain);
% g = merge(g, idx.'); 

% [TODO]: Do we want to do this?
g = simplify(g);

end
