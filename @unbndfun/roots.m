function r = roots(f, varargin)
%ROOTS   Roots of an UNBNDFUN in an unbounded domain.
%   ROOTS(F) returns the real roots of the UNBNDFUN F.
%
%   ROOTS(F, OPTIONS) modifies the default ROOTS properties, by passing the
%   OPTIONS to the rootfinding method of the ONEFUN of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call ROOTS@CLASSICFUN:
r = roots@classicfun(f, varargin{:});

% Try to get rid of spurious roots which are caused by fast decay of the
% function at +/- Inf.

ends = get(f, 'domain');
numRoots = length(r);

% If f is defined on an interval extending to -Inf and, given a root r(n), it
% is "essentially zero" on the subinterval [-Inf, r(n)], the root is probably
% spurious, so mark it for filtering.
if ( isinf(ends(1)) )
    for n = 1:1:numRoots
        rf = restrict(f, [-Inf, r(n)]);
        if (normest(rf) < get(f, 'vscale').*get(f, 'epslevel'))
            mask(n) = true;
        else
            break;
        end
    end
end

% Same thing, but for functions defined on an interval extending to +Inf.
if ( isinf(ends(2)) )
    for n = numRoots:-1:1
        rf = restrict(f, [r(n), Inf]);
        if (normest(rf) < get(f, 'vscale').*get(f, 'epslevel'))
            mask(numRoots - n + 1) = true;
        else
            break;
        end
    end
end

r(mask) = [];

end
