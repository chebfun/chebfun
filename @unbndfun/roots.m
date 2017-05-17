function r = roots(f, varargin)
%ROOTS   Roots of an UNBNDFUN in an unbounded domain.
%   ROOTS(F) returns the real roots of the UNBNDFUN F.
%
%   ROOTS(F, OPTIONS) modifies the default ROOTS properties, by passing the
%   OPTIONS to the rootfinding method of the ONEFUN of F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call ROOTS@FUN.  Roots of an UNBNDFUN should always be pruned to try and
% remove spurious roots caused by fast decay of f near infinite endpoints.
if ( (nargin > 1) && isstruct(varargin{1}) )
    varargin{1}.filter = @unbndFilter;
    r = roots@classicfun(f, varargin{:});
else
    r = roots@classicfun(f, varargin{:}, 'filter', @unbndFilter);
end

end

function r = unbndFilter(r, f)
%UNBNDFILTER   Try to detect and remove spurious roots near infinity.
%   This function exists primarily to help with filtering of spurious roots for
%   functions represented using UNBNDFUN which decay at +/-Inf. 
%
%   The idea is that if a function is near zero for all values between a root
%   and the nearest interval endpoint, its more likely that the root which
%   emerged was spurious, showing up only due to rounding error in a function
%   which decays towards +/-Inf.
%
%   Since this filtering is actually performed at the TECH level, we look at the
%   values of the function at some points in [-1, r] and [r, 1].

% TODO: This assumes that the tech eventually used to find the roots lives on
% [-1,1] and has a .vscale property. Currently this is only a verbal agreement,
% as there is no abstract tech class.

numRoots = length(r);
mask = false(numRoots, 1);
tol = 10*eps*get(f, 'vscale');

% We sample at an arbitrary number of points between the located root and
% the nearest enpoint. We take max(20, numRoots), with the reasoning being that
% if there are many roots, then the function has some complex oscillatory
% behaviour that we might need to capture.
n = max(numRoots, 20);

% Filter the roots at the left endpoint.
if ( abs(feval(f, -1)) < tol )
    for k = 1:1:numRoots
        if ( r(k) > 0 )
            % This roots is far from -1.
            continue
        end
        % Equispaced grid:
        testGrid = linspace(-1, r(k), n);
        if ( norm(feval(f, testGrid), Inf) < tol )
            % Remove this root.
            mask(k) = true;
        else
            % We needn't check subsequent roots.
            break
        end
    end
end

% Filter the roots at the right endpoint.
if ( abs(feval(f, 1)) < tol )
    for k = numRoots:-1:1
        if ( r(k) < 0 )
            % This roots is far from +1.
            continue
        end
        % Equispaced grid:
        testGrid = linspace(r(k), 1, n);
        if ( norm(feval(f, testGrid), Inf) < tol )
            % Remove this root.
            mask(k) = true;
        else
            % We needn't check subsequent roots.
            break
        end
    end
end

r(mask) = [];

end
