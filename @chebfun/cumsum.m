function f = cumsum(f, m, dim)
%CUMSUM   Indefinite integral of a CHEBFUN.
%   G = CUMSUM(F) is the indefinite integral of the column CHEBFUN F. G will
%   typically be normalised so that G(F.domain(1)) = 0.
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

% [TODO]: Update the above help text once we have singfun.
%   G will typically be normalised so that G(F.domain(1)) = 0.  The exception to
%   this is when computing indefinite integrals of functions which are not
%   integrable at the left boundary. In this case, the arbitrary constant in the
%   indefinite integral is chosen to make the representation of G as simple as
%   possible. Dirac deltas already existing in F will decrease their degree.

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

% Loop m times:
for l = 1:m
    
    % Get the level 2 (delta function) impulse data:
    if ( size(f.impulses, 3) > 1 )
        deltas = f.impulses(:,:,2);
    else
        deltas = zeros(length(dom), numCols);
    end
    
    fa = deltas(1,:);
    for j = 1:numFuns
        
        cumsumFunJ = cumsum(funs{j});

        if ( numFuns > 1 )
            % We evaluate cumsumFunJ at both endpoints to see if the function
            % value is finite. Infinite function value at the endpoints means
            % cumsumFunJ is a singfun. If the function is infinite at the left
            % end, then it doesn't make sense to shift the function up or down 
            % to make its value zero at the left endpoint. If the function is 
            % infinite at the right end, then shifting function up or down means
            % we are computing SINGFUN + constant, which, in general, is not 
            % accurate running out of the points which are allowed. Such a 
            % difficulty may be mitigated when SING MAP is re-adopted. Also note
            % that we can't prevent the case where a constant is added to a
            % SINGFUN with positive exponent, e.g. x.^(0.5) + 1, since in this
            % level we don't know the type of the ONEFUN of each FUN.
            
            lval = get(cumsumFunJ, 'lval');
            rval = get(cumsumFunJ, 'rval');
            if ( all(isfinite(lval)) && all(isfinite(rval)) )
                cumsumFunJ = cumsumFunJ - lval;
                if ( all(fa ~= 0) && all(isfinite(fa)) )
                    cumsumFunJ = cumsumFunJ + fa;
                end
                fa = get(cumsumFunJ, 'rval') + deltas(j+1,:); 
            end
        end
        
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

