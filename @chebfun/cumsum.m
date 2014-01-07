function f = cumsum(f, m, dim)
%CUMSUM   Indefinite integral of a CHEBFUN.
%   G = CUMSUM(F) is the indefinite integral of the column CHEBFUN F. G will
%   typically be normalised so that G(F.domain(1)) = 0. The exception to this is
%   when computing indefinite integrals of functions whose indefinite integrals
%   have singularities. In this case, the arbitrary constant in the indefinite
%   integral is chosen to make the representation of G as simple as possible.
%
%   CUMSUM(F, N) returns the Nth integral of F. If N is not an integer then
%   CUMSUM(F, N) returns the fractional integral of order N as defined by the
%   Riemann-Liouville integral.
%
%   CUMSUM(F, N, 2) will take the Mth cumulative sum over the columns F an
%   array-valued BNDFUN.
%
% See also SUM, INTEGRAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Update the above help text once we have deltafun. Dirac deltas already
% existing in F will decrease their degree.

% Trivial case:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    m = 1;
end
if ( nargin < 3 )
    % Continuous dimension by default:
    dim = 1 + f.isTransposed;
end

if ( round(m) ~= m )
    % Fractional integral:
    % [TODO]: Implement this!
    error('CHEBFUN:cumsum:notImplemented', ...
        'Fractional antiderivatives not yet implemented.');
    f = fracCalc(f, m);
    return
end

if ( ( dim == 1 && ~f.isTransposed ) || ( dim == 2 && f.isTransposed ) )
    % Continuous dimension:
    f = cumsumContinousDim(f, m);
else
    % Finite dimension:
    f = cumsumFiniteDim(f, m);
end

end

function f = cumsumContinousDim(f, m)
% CUMSUM over continuous dimension.

% Get some basic information from f:
dom = f.domain;
funs = f.funs;
numFuns = numel(funs);
numCols = size(f.funs{1}, 2);

% Preprocess for singular cases. If a FUN has nontrivial exponents at both
% endpoints, then a break point is introduced to accommodate the disability
% of @SINGFUN/CUMSUM for handling functions with non-zero exponent at both
% ends.
domOld = f.domain;
breakPoints = [];
toBreak = 0;

% Is there any SINGFUN involved has non-zero exponent at both ends? If so,
% set the flag for introducing new break points true and then take the
% mid-point of these domains as the new break points.
for j = 1:numFuns
    if ( isa(f.funs{j}.onefun, 'singfun') && all(f.funs{j}.onefun.exponents) )
        toBreak = 1;
        breakPoints = [breakPoints mean(f.domain(j:j+1))];
    end
end

% Introduce new break points using RESTRICT.
if ( toBreak )
    domNew = sort([domOld, breakPoints]);
    f = restrict(f, domNew);
    dom = f.domain;
    funs = f.funs;
    numFuns = numel(funs);
end

% Loop m times:
for l = 1:m
    
    % Get the level 2 (delta function) impulse data:
    if ( size(f.impulses, 3) > 1 )
        deltas = f.impulses(:,:,2);
    else
        deltas = zeros(length(dom), numCols);
    end
    
    rval = deltas(1,:);
    
    % Main loop for looping over each piece and do the integration:
    for j = 1:numFuns
        
        % In fact, CUMSUM@BNDFUN will check if the current piece, i.e.
        % cumsumFunJ.onefun is a SINGFUN. If so, then we don't want to shift the
        % current piece up or down to stick the left end of the current piece to the
        % right end of the last one, since SINGFUN + CONSTANT won't be accurate and
        % may trigger annoying SINGFUN warning messages. Such a difficulty may
        % be mitigated when SING MAP is re-adopted. Also if the last piece is
        % infinite at the right end, then shifting the current piece to concatenate
        % doesn't make any sense.
        
        % Call CUMSUM@BNDFUN:
        cumsumFunJ = cumsum(funs{j}, 1, 1, rval);
        
        % Update the value of the right end:
        rval = get(cumsumFunJ, 'rval') + deltas(j+1,:);
        
        % Store the current piece:
        funs{j} = cumsumFunJ;
        
    end
    
    % Get the new impulse data:
    newImps = chebfun.getValuesAtBreakpoints(funs, dom);
    f.impulses = cat(3, newImps, f.impulses(:,:,3:end));
    
end

% Append the updated FUNs:
f.funs = funs;

end

function f = cumsumFiniteDim(f, m)
% CUMSUM over finite dimension.

for k = 1:numel(f.funs)
    f.funs{k} = cumsum(f.funs{k}, m, 2);
end

end
