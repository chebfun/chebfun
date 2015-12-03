function out = isequal(f, g)
%ISEQUAL   Equality test for two CHEBFUNs.
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the CHEBFUN objects F and G
%   contain identical breakpoints and FUNS, and logical 0 (FALSE) otherwise.
%
% See also EQ.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check the empty case:
if ( isempty(f) && isempty(g) )
    out = true;
    return
end

% Check the domains, transpose states, and the number of columns match:
if ( ~domainCheck(f, g) || f(1).isTransposed ~= g(1).isTransposed || ...
    numColumns(f) ~= numColumns(g) )
    out = false;
    return
end

if ( numel(f) == 1 && numel(g) == 1 )
    % Array-valued CHEBFUN case:
    out = columnIsequal(f, g);
else
    % QUASIMATRIX case:
    f = mat2cell(f);
    g = mat2cell(g);
    % Loop over the columns:
    for k = 1:numel(f)
        out = columnIsequal(f{k}, g{k});
        if ( ~out ), break, end
    end
end

end

function out = columnIsequal(f, g)

% Check the number of FUNs:
if ( numel(f.funs) ~= numel(g.funs) )
    out = false;
    return
end

% Check the pointValues:
tol = 1e1*eps;
if ( norm(f.pointValues(:) - g.pointValues(:), inf) > tol )
    out = false;
    return
end

% Check the FUNs:
for j = 1:numel(f.funs)
    if ( ~isequal(f.funs{j}, g.funs{j}) )
        out = false;
        return
    end
end

% If everything matched up, f and g must be the same:
out = true;

end
