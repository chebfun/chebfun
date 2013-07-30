function out = isequal(f, g)
%ISEQUAL    Equality test for two chebfuns.
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the CHEBFUN objects F and G
%   contain identical breakpoints and funs, and logical 0 (FALSE) otherwise.
%
% See also chebfun/eq.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check the domains match:
if ( numel(f.domain) ~= numel(g.domain) || ...
        norm(f.domain - g.domain, inf) > hscale(f)*eps );
    out = false;
    return
end

% Check the impulses:
if ( norm(f.impulses - g.impulses, inf) > get(f, 'epslevel') )
    out = false;
    return
end

% Check the funs:
for j = 1:numel(f.funs)
    if ( ~isequal(f.funs{j}, g.funs{j}) )
        out = false;
        return
    end
end

% If everything matched up, f and g must be the same:
out = true;

end