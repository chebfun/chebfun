function f = diff(f, n, dim)
%DIFF   Differentiation of a CHEBFUN.
%   DIFF(F) is the derivative of the columnn chebfun F. At discontinuities, DIFF
%   creates a Dirac delta with coefficient equal to the size of the jump. Dirac
%   deltas already existing in F will increase their degree. 
%
%   DIFF(F, N) is the Nth derivative of F.
%
%   DIFF(F, N, 2) of an array-valued row CHEBFUN (or DIFF(F, N, 1) of a
%   array-valued column CHEBFUN) computes the Nth difference accross the columns
%   (or rows) of F.
%
% See also SUM, CUMSUM.

%   [TODO]: Fractional derivatives. DIFF(F, ALPHA), when ALPHA is not an
%   integer, offers some support for fractional derivatives (of degree ALPHA) of
%   F. For ALPHA > 1 the Riemann- Liouville definition is used by default. On
%   can switch to the Caputo definition with a call of the form DIFF(F, ALPHA,
%   'Caputo'). [Requires SINGFUN].

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Trivial case:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    n = 1;
end
if ( nargin < 3 )
    dim = 1;
end

if ( ~any(dim == [1, 2]))
    error('CHEBFUN:cumsum:dim', 'Index exceeds matrix dimensions.');
elseif ( round(dim) ~= dim )
    error('CHEBFUN:cumsum:dim', 'Dimension must either be 1 or 2.');
end
    
if ( round(n) ~= n )
    % Fractional integral:
    % [TODO]: Implement this!
    f = fracCalc(f, n);
    return
end

% Diff across columns (or rows for a transposed) array-valued chebfun:
if ( xor(f.isTransposed, dim == 2) )
    for k = 1:numel(f.funs)
        f.funs{k} = diff(f.funs{k}, n, 2);
    end
    return
end

% Grab some fields from f:
funs = f.funs;
numFuns = numel(funs);
numCols = size(f.funs{1}, 2);
imps = f.impulses;

% Set a tolerance: (used for introducing Dirac deltas at jumps)
tol = epslevel(f)*hscale(f);

% Loop n times for nth derivative:
for j = 1:n
    vs = get(f, 'vscale-local');

    % Detect jumps in the original function and create new deltas.
    newDeltas = zeros(numFuns + 1, numCols, 1);
    for k = 1:numFuns-1
        jmp = get(funs{k+1}, 'lval') - get(funs{k}, 'rval');
        scl = 0.5*(vs(k) + vs(k+1,:));
        if ( any( abs(jmp) > tol*scl ) )
           newDeltas(k+1,:,:) = jmp;
        end
    end

    % Differentiate each FUN in turn:
    for k = 1:numFuns
        funs{k} = diff(funs{k});
    end

    % Compute new function values at breaks using JUMPVALS():
    imps(:,:,1) = chebfun.getValuesAtBreakpoints(funs);
    % Update impulses:
    if ( size(imps, 3) > 1 )
       imps = cat(3, imps(:,:,1), newDeltas, imps(:,:,2:end));
    elseif ( any(newDeltas) )
       imps(:,:,2) = newDeltas;
    end

    % Reassign data to f:
    f.funs = funs;
    f.impulses = imps;

end

end