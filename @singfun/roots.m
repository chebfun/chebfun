function out = roots(f, varargin)
%ROOTS   Roots of a SINGFUN in the interval [-1,1].
%   ROOTS(F) returns the real roots of the SINGFUN F in the interval [-1,1].
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs (strings) and VALs may be any of the following:
%
%   ALL: 
%       [0] - Return only real-valued roots in [-1,1].
%        1  - Return roots outside of [-1,1] (including complex roots).
%
%   RECURSE:
%        0  - Compute roots without interval subdivision (slower).
%       [1] - Subdivide until length(F) < 50. (causes additional complex roots).
%
%   PRUNE:
%       [0]
%        1  - Prune 'spurious' complex roots if ALL == 1 and RECURSE == 0.
%
% Note that the roots of the smoothPart is computed by calling ROOTS in
% SMOOTHFUN. See ROOTS in SMOOTHFUN and levels below for more documentations.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Set up a tolerance:
tol = eps;

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

% Compute roots of the smooth part:
out = roots(f.smoothPart, varargin{:} );
anyRoots = ~isempty(out);

% Determine if the endpoints are roots:
if ( f.exponents(1) > 0 )
    if ( anyRoots )
        if ( abs(1 + out(1)) < tol )
            out(1) = -1;
        else
            out = [-1; out];
        end
    else
        out = -1;
    end
end

if ( f.exponents(2) > 0 )
    if ( anyRoots )
        if ( abs(1 - out(end)) < tol )
            out(end) = 1;
        else
            out = [out; 1];
        end
    else
        out = [out; 1];
    end
end

end
