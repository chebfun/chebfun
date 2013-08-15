function out = isequal(f, g)
%ISEQUAL   Equality test for two CHEBFUNs.
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the CHEBFUN objects F and G
%   contain identical breakpoints and funs, and logical 0 (FALSE) otherwise.
%
% See also EQ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check the empty case:
if ( isempty(f) && isempty(g) )
    out = true;
    return
end

% Otherwise, intialise output as false:
out = false;

% Check the domains, sizes, and transpose states match:
if ( ~domainCheck(f, g) || ~isequal(size(f), size(g)) || ...
     ~isequal(f.isTransposed, g.isTransposed) )
    return
end

% Check the impulses:
tol = max(epslevel(f), epslevel(g));
if ( norm(f.impulses - g.impulses, inf) > tol )
    return
end

% Check the FUNs:
for j = 1:numel(f.funs)
    if ( ~isequal(f.funs{j}, g.funs{j}) )
        return
    end
end

% If everything matched up, f and g must be the same:
out = true;

end
