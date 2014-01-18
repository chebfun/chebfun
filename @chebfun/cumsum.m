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
%   array-valued CHEBFUN or quasimatrix.
%
% See also SUM, INTEGRAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Update the above help text once we have deltafun. Dirac deltas already
% existing in F will decrease their degree.

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
    % [TODO]: Implement this!
    error('CHEBFUN:cumsum:notImplemented', ...
        'Fractional antiderivatives not yet implemented.');
    f = fracCalc(f, m);
    return
end

if ( ( dim == 1 && ~f(1).isTransposed ) || ( dim == 2 && f(1).isTransposed ) )
    % Continuous dimension:
    for k = 1:numel(f)
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
numCols = size(f.funs{1}, 2);
dom = f.domain;
numFuns = numel(f.funs);

% Loop m times:
for l = 1:m
    
    % Get the level 2 (delta function) impulse data:
    if ( size(f.impulses, 3) > 1 )
        deltas = f.impulses(:,:,2);
    else
        deltas = zeros(length(dom), numCols);
    end
    
    rval = deltas(1,:);
    funs = [];
    
    % Main loop for looping over each piece and do the integration:
    for j = 1:numFuns
        
        % CUMSUM@BNDFUN will check if the current piece, i.e. cumsumFunJ.onefun 
        % is a SINGFUN. If so, then we don't want to shift the current piece up 
        % or down to stick the left end of the current piece to the right end of
        % the last one, since SINGFUN + CONSTANT won't be accurate and may 
        % trigger annoying SINGFUN warning messages. Such a difficulty may be 
        % mitigated when SING MAP is re-adopted. Also if the last piece is
        % infinite at the right end, then shifting the current piece to 
        % concatenate doesn't make any sense.
        
        % Call CUMSUM@BNDFUN:
        cumsumFunJ = cumsum(f.funs{j}, 1, 1, rval);
        
        % [TODO]: Check why deltas appears here. 
        
        if ( iscell( cumsumFunJ ) )
            % Update the value of the right end:
            rval = get(cumsumFunJ{2}, 'rval') + deltas(j+1,:);
        else
            % Update the value of the right end:
            rval = get(cumsumFunJ, 'rval') + deltas(j+1,:);
        end
        
        % Store FUNs:
        funs = [funs, {cumsumFunJ}];
        
    end
    
    % Get the new impulse data:
    newImps = chebfun.getValuesAtBreakpoints(funs, dom);
    f.impulses = cat(3, newImps, f.impulses(:,:,3:end));
    
    % Append the updated FUNs:
    f.funs = funs;
    
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
