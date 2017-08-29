function h = eq(f, g)
%==   Equality operator for CHEBFUN objects.
%   H = F == G, where F and/or G are CHEBFUN objects, constructs a logical
%   CHEBFUN H which is true (i.e., takes the value 1) where F == G, and false
%   (0) elsewhere.
%
% See also NE, ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case:
if ( isempty(f) )
    h = f;
    return
elseif ( isempty(g) )
    h = g;
    return
end

% Array-valued?
if ( min(size(f)) > 1 || min(size(g)) > 1 )
	error('CHEBFUN:CHEBFUN:eq:array', ...
        '== does not support array-valued CHEBFUN objects.');
end

% Call SIGN() to do the work:
h = sign(f - g);

% Get value in interior of each FUN by taking left-sided limit at the breaks:
vals = get(h, 'rval-local');

% Set FUNs that are 0 to 1:
for k = 1:numel(h.funs)
    if ( vals(k) == 0 )
        h.funs{k} = 1 + 0*h.funs{k};
    else
        h.funs{k} = 0*h.funs{k};
    end
end

% pointValues:
nanMask = isnan(h.pointValues);
ind = ~nanMask;
h.pointValues(ind) = ~h.pointValues(ind);

% Tidy the result:
h = merge(h);

end
