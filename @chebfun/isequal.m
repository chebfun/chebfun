function out = isequal(f, g)
%ISEQUAL   Equality test for two CHEBFUNs.
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the CHEBFUN objects F and G
%   contain identical breakpoints and funs, and logical 0 (FALSE) otherwise.
%
% See also CHEBFUN/EQ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Intialise output as false:
out = false;

% Check the domains match:
if ( ~domainCheck(f, g) );
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
