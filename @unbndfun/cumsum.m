function [f, rVal] = cumsum(f, k, dim)
%CUMSUM   Indefinite integral of an UNBNDFUN.
%   CUMSUM(F) is the indefinite integral of the UNBNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, K) will compute the Kth indefinite integral with the constant of
%   integration chosen so that each intermediate integral evaluates to 0 at x=a.
%   Thus CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, K, 2) will take the Kth cumulative sum over the columns F an
%   array-valued UNBNDFUN.
%
%   [F, RVAL] = CUMSUM(F), [F, RVAL] = CUMSUM(F, K), and [F, RVAL] = 
%   CUMSUM(F, K, 2) will do the same thing as above, but also return the value
%   of the integral at the right endpoint, which will be used at CHEBFUN level
%   for concatenating neighboring pieces.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trivial case of an empty UNBNDFUN:
if ( isempty(f) )
    return
end

%% Parse inputs:
if ( nargin == 1 || isempty(k) )
    % Compute first indefinite intergral by default.
    k = 1;
end

if ( nargin < 3 )
    % Assume dim = 1 by default.
    dim = 1;
end

%% Treat separately for different dimensions:

if ( dim == 1 )
    
    % Compute the indefinite integral along the continuous dimension:
    for j = 1:k
        if ( iscell(f) )
            % If F turns out to be a cell of two pieces, we integrate each
            % piece separately:
            f{1} = cumsumCtsDim(f{1});
            f{2} = cumsumCtsDim(f{2});
        else
            f = cumsumCtsDim(f);
        end
    end
    
elseif ( dim == 2 )
    
    % When dim = 2, we compute the cumlative sum over columns:
    f.onefun = cumsum(f.onefun, k, 2);
    
else
    error('CHEBFUN:UNBNDFUN:cumsum:input', ...
        'The third argument is unrecognizable.');
end

% Value of F at the right endpoint:
if ( iscell(f) )
    rVal = cellfun(@(f) get(f, 'rval'), f);
else
    rVal = get(f, 'rval');
end

end

function g = cumsumCtsDim(f, pref)

% Make a copy of F:
g = f;

% Rescaling factor is the derivative of the forward map:
pref.singPrefs.exponents = g.mapping.forDerExps;
rescaleFactor = onefun.constructor(@(x) g.mapping.forDer(x), [], [], pref);
numRoots = -repmat(pref.singPrefs.exponents.', 1, size(g, 2));

% Try to see if we can extract boundary roots:
[h, rootsLeft, rootsRight] = extractBoundaryRoots(g.onefun, numRoots);

if ( all(rootsLeft == numRoots(1,:)) && all(rootsRight == numRoots(2,:)) )
    
    % The ONEFUN of the integral of F should be the integral of the ONEFUN
    % of the F multiplied by the derivative of the forward map. Here the
    % singularities of the RESCALEFACTOR is cancelled off by the boundary
    % roots of H. Therefore, only the smoothPart of RESCALEFACTOR is
    % involved.
    g.onefun = cumsum(h.*rescaleFactor.smoothPart);
    
else
    
    % The ONEFUN of the integral of F should be the integral of the
    % ONEFUN of the F multiplied by the derivative of the forward
    % map.
    g.onefun = g.onefun.*rescaleFactor;
    g.onefun = cumsum(g.onefun);
    
end

end