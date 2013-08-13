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
%   If F is a array-valued CHEBTECH then there is no reason to expect each
%   column to have the same number of roots. In order to return a useful output,
%   the roots of each column are computed and then padded with NaNs so that a
%   matrix may be returned. The columns of R = ROOTS(F) correspond to the
%   columns of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Set up a tolerance:
tol = get(f.smoothPart, 'epslevel')*get(f.smoothPart, 'vscale');

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

out = roots(f.smoothPart, varargin{:} );
anyRoots = ~isempty(out);

% Grap the exponents:
exps = f.exponents;

% Determine if the endpoints are roots:
if ( exps(1) > 0 )
    if ( anyRoots )
        if ( abs(1+out(1)) < tol )
            out(1) = -1;
        else
            out = [-1; out];
        end
    else
        out = -1;
    end
end

if ( exps(2) > 0 )
    if ( anyRoots )
        if ( abs(1-out(end)) < tol )
            out(end) = 1;
        else
            out = [out; 1];
        end
    else
        out = [out; 1];
    end
end

end