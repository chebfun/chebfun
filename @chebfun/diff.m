function F = diff(F, n, dim)
%DIFF   Differentiation of a CHEBFUN.
%   DIFF(F), when F is a column CHEBFUN, computes a column CHEBFUN whose columns
%   are the derivatives of the corresponding columns in F.  At discontinuities,
%   DIFF creates a Dirac delta with coefficient equal to the size of the jump.
%   Dirac deltas already existing in F will increase their degree.
%
%   DIFF(F), when F is an array-valued row CHEBFUN or a quasimatrix, computes
%   the first-order finite difference of F along its rows. The resulting row
%   CHEBFUN will have one row less than the number of rows in F.
%
%   DIFF(F, N) or DIFF(F, N, 1) computes the Nth derivative of F if F is a
%   column CHEBFUN and the Nth-order finite difference of F along its rows if F
%   is a row CHEBFUN.
%
%   DIFF(F, N, 2) is the Nth-order finite difference of F along its columns if
%   F is a column CHEBFUN and the Nth derivative of F if F is a row CHEBFUN.
%
%   DIFF(F, MU), when MU is not an integer returns the MUth Riemann-Liouville
%   fractional derivative of the CHEBFUN F. DIFF(F, MU, 'Caputo') uses instead
%   the Caputo definition. See [1] for definitions. In either case, an error is
%   thrown if F is not smooth or is defined on an unbounded domain.
%
% See also SUM, CUMSUM.

% References:
%  [1] http://en.wikipedia.org/wiki/Fractional_calculus#Fractional_derivatives

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( isempty(F) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    n = 1;
end
if ( nargin < 3 )
    dim = 1;
end

if ( isnumeric(dim) && ~any(dim == [1, 2]) )
    error('CHEBFUN:CHEBFUN:diff:dim', 'Dimension must either be 1 or 2.');
end
    
if ( round(n) ~= n )
    % Fractional derivative:
    F = fracDiff(F, n, dim);
    return
end

if ( xor(F(1).isTransposed, dim == 2) )
    % Diff across columns (or rows for a transposed) array-valued CHEBFUN:
    F = diffFiniteDim(F, n);
else
    % Diff along continuous dimension (i.e., dF/dx):
    for k = 1:numel(F)
        F(k) = diffContinuousDim(F(k), n);
    end
end

end

function f = diffFiniteDim(f, n)
% Differentiate across the finite dimension (i.e., across columns).
if ( numel(f) == 1 )
    % Array-valued CHEBFUN case:
    for k = 1:numel(f.funs)
        f.funs{k} = diff(f.funs{k}, n, 2);
    end
else
    % Quasimatrix case:
    numCols = numel(f);
    if ( numCols <= n )
        f = chebfun();
    else
        for j = 1:n
            for k = 1:numCols-j
                f(k) = f(k+1) - f(k);
            end
        end
    end
    f = f(1:numCols-n);
end

end

function f = diffContinuousDim(f, n)
% Differentiate along continuous dimension (i.e., df/dx).

% Grab some fields from f:
funs = f.funs;
numFuns = numel(funs);
numCols = size(f.funs{1}, 2);

% Grb a preference:
pref = chebfunpref();

% Array for keeping track of point values in case there is a delta function at
% a break point.
infVals = [];

if ( ~pref.enableDeltaFunctions )
    % No delta functions:
    
    % Differentiate each FUN in turn:
    for k = 1:numFuns
        funs{k} = diff(funs{k}, n);
    end
  
else
    % Delta functions are enabled:
    
    % Tolerance used for introducing Dirac deltas at jumps:
    deltaTol = pref.deltaPrefs.deltaTol;
    
    % Loop n times for nth derivative:
    for j = 1:n    
        % Detect jumps in the original function and create new deltas.
        deltaMag = getDeltaMag();
        infVals = inf * deltaMag;
        % Differentiate each FUN in turn:
        for k = 1:numFuns
            funs{k} = diff(funs{k});
            % If there is a delta function at the join, recreate the FUN using the
            % DELTAFUN constructor:
            funs{k} = makeDeltaFun(funs{k}, deltaMag(k:k+1,:));
        end

    end         
    
end

% Compute new function values at breaks:
pointValues = chebfun.getValuesAtBreakpoints(funs);
I = isinf(infVals);
pointValues(I) = infVals(I); 

% Reassign data to f:
f.funs = funs;
f.pointValues = pointValues;

    function deltaMag = getDeltaMag()
        deltaMag = zeros(numFuns + 1, numCols);
        % Loop through the funs:
        for l = 1:(numFuns - 1)
            % Extract the jump vlaues between two funs:
            jmp = get(funs{l+1}, 'lval') - get(funs{l}, 'rval');
            % Assign these jumps to the deltaMag matrix:
            if ( any(abs(jmp) > deltaTol ) )
                deltaMag(l+1, :) = jmp;
            end
        end
    end

    function f = makeDeltaFun(f, deltaMag)
        if ( any(abs(deltaMag(:)) > deltaTol) )
            % [TODO]: This does not handle array-valuedness at the moment.
            if ( size(deltaMag, 2) > 1 )
                warning('CHEBFUN:CHEBFUN:diff:diffContinuousDim:makeDeltaFun:array', 'No support here for array-valuedness at the moment.');
                deltaMag = [0 ; 0];
            end
            % New delta functions are only possible at the ends of the domain:
            data.domain = f.domain;
            data.deltaMag = deltaMag.'/2;
            data.deltaLoc = f.domain;
            % Add new delta functions to the existing fun:
            f = deltafun(0, data, pref) + f;
        end
    end

end
