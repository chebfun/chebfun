function h = le(f, g)
% LE (<=) for CHEBFUN object.
%   H = F <= G, where F and/or G are CHEBFUN objects, constructs a logical
%   CHEBFUN H which is true (i.e., takes the value 1) where F <= G, and false
%   (0) elsewhere.
%
% See also GE, LT, GT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
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
	error('CHEBFUN:le:array', ...
        '<= does not suport array-valued CHEBFUN objects.');
end

% Call HEAVISIDE() to do the work:
h = heaviside(g - f);

% Get a value in the interior of each FUN:
% % Evaluate at middle of domain:
% x = (h.domain(1:end-1) + h.domain(2:end)).'/2;
% vals = feval(h, x);
% Take the left-sided limit at the breaks:
vals = get(h, 'rval-local');

% Set FUNs that are 0.5 to 0 FUNs:
for k = 1:numel(h.funs)
    if ( vals(k) == .5 )
        h.funs{k} = 1 + 0*h.funs{k};
    end
end

% Set impulses that are 0.5 to 0 FUNs:
h.impulses(h.impulses == .5) = 1;

% Tidy the result:
h = merge(h);

end
