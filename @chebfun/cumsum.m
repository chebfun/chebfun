function f = cumsum(f, m, dim)
%CUMSUM   Indefinite integral of a CHEBFUN.
%   G = CUMSUM(F) is the indefinite integral of the column CHEBFUN F. G will
%   typically be normalised so that G(F.domain(1)) = 0. The exception to this is
%   when computing indefinite integrals of functions whose indefinite integrals
%   have singularities. In this case, the arbitrary constant in the indefinite
%   integral is chosen to make the representation of G as simple as possible.
%
%   CUMSUM(F, N) returns the Nth integral of F. If N is not an integer then
%   CUMSUM(F, MU) returns the fractional integral of order MU as defined by the
%   Riemann-Liouville integral [1]. In the latter case an error is thrown if F
%   is not smooth or is defined on an unbounded domain.
%
%   CUMSUM(F, N, 2) will take the Nth cumulative sum over the columns F an
%   array-valued CHEBFUN or quasimatrix.
%
% See also SUM, INTEGRAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% References:
%  [1] http://en.wikipedia.org/wiki/Riemann%E2%80%93Liouville_integral

% TODO: The input sequence is not the same as MATLAB. In particular, MATLAB only
% supports m = 1.

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
    dim = 1 + f(1).isTransposed;
end

if ( round(m) ~= m )
    % Fractional integral:
    f = fracInt(f, m);
    return
end

if ( ( dim == 1 && ~f(1).isTransposed ) || ( dim == 2 && f(1).isTransposed ) )
    % Continuous dimension:
    for k = 1:numel(f)
        % Handle delta functions:
        if ( isdelta(f(k)) )
            f(k) = pruneDeltas(f(k));
        end
        f(k) = cumsumContinousDim(f(k), m);
    end
else
    % Finite dimension:
    f = cumsumFiniteDim(f, m);
end

end

function f = cumsumContinousDim(f, m)
% CUMSUM over continuous dimension.

% Get some basic information from f:
numFuns = numel(f.funs);

transState = f(1).isTransposed;

% If f is periodic, i.e. based on a periodic TECH check its mean:
if ( isPeriodicTech(f) )
    % Obtain Fourier coefficients {c_k}
    c = trigcoeffs(f); 
    numCoeffs = size(c, 1);
    % index of constant coefficient
    if ( mod(numCoeffs, 2) == 0 )
       ind = numCoeffs/2 + 1;
    else
       ind = (numCoeffs + 1)/2;
    end
    if ( any(abs(c(ind,:)) > 1e1*vscale(f)*eps) )
        % Mean is not zero, convert it to a CHEBTECH based chebfun:
        f = chebfun(f);
    end
end

% Loop m times:
for l = 1:m

    funs = {};
    rVal = 0;   
    % Main loop for looping over each piece and do the integration:
    for j = 1:numFuns

        % Call FUN/CUMSUM():        
        newFuns = cumsum(f.funs{j});
               
        if ( ~iscell( newFuns ) )
            newFuns = {newFuns};
        end
                        
        % Add the constant term that came from the left:
        for k = 1:numel(newFuns)
            newFuns{k} = newFuns{k} + rVal;
        end

        % Store FUNs:
        funs = [funs, newFuns]; %#ok<AGROW>
        
                    
        % Update the rval:
        rVal = get(funs{end}, 'rval');
        
    end
    
    % Assemble the new CHEBFUN:
    f = chebfun(funs);
    
end

if ( transState )
    f = f.';
end

end

function f = cumsumFiniteDim(f, m)
% CUMSUM over finite dimension.

if ( numel(f) == 1 )
    for k = 1:numel(f.funs)
        f.funs{k} = cumsum(f.funs{k}, m, 2);
    end
else
    numCols = numel(f);
    for j = 1:m
        for k = 2:numCols-j
            f(k) = f(k) + f(k-1);
        end
    end
end

end

function f = pruneDeltas(f)
%PRUNEDELTAS   Transfers deltas at the right endpoint of funs to the next fun.
%   F is assumed to be ba chebfun with a single column.

% Go through each fun in pairs:
nFuns = numel(f.funs);
for k = 1:(nFuns-1);
    if ( isdelta(f.funs{k}) )
        [f.funs{k}, f.funs{k+1}] = transferDeltas(f.funs{k}, f.funs{k+1});
    end
end
end
