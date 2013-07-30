function f = diff(f, n, pref)
%DIFF   Differentiation of a chebfun.
%   DIFF(F) is the derivative of the columnn chebfun F. At discontinuities, DIFF
%   creates a Dirac delta with coefficient equal to the size of the jump. Dirac
%   deltas already existing in F will increase their degree. DIFF(F, N) is the
%   Nth derivative of F.
%
%   DIFF(F, N, 2) of an array-valued row chebfun (or DIFF(F, N, 1) of a
%   array-valued column chebfun) computes the Nth difference accross the columns
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
    pref = chebfun.pref;
elseif ( nargin == 2 )
    if ( isstruct(n) )
        pref = n;
        n = 1;
    else
        pref = chebfun.pref;
    end
end

% Diff across columns (or rows for a transposed) array-valued chebfun:
if ( xor(f.isTransposed, ~isstruct(pref) && pref == 2) )
    for k = 1:numel(f.funs)
        f.funs{k} = diff(f.funs{k}, n, 2);
    end
    return
end

% Set a tolerance: (used for introducing Dirac deltas at jumps)
tol = max(10*pref.chebfun.eps, 1e-12);

% Grab some fields from f:
funs = f.funs;
nfuns = numel(funs);
imps = f.impulses;

if ( numel(funs) == 1 )
    % Single piece (no new impulses).

    % Can simply call FUN/DIFF() for nth derivative:
    funs{1} = diff(funs{1}, n);
    % Update the impulses:
    imps = [zeros(n,1) ; imps];

else
    % Multiple funs.

    % Loop n times for nth derivative:
    for j = 1:n

        % Detect jumps in the original function and create new impulses.
        newImps = zeros(1, nfuns + 1, size(funs{1}, 2));
        for k = 1:nfuns-1
            jmp = get(funs{k+1}, 'lval') - get(funs{k}, 'rval');
            scl = .5*(get(funs{k}, 'vscale') + get(funs{k+1}, 'vscale'));
            if any( abs(jmp) > tol*scl )
               newImps(:,k+1,:) = jmp;
            end
        end

        % Differentiate each fun in turn:
        for k = 1:nfuns
            funs{k} = diff(funs{k});
        end

        % Compute new function values at breaks using JUMPVALS():
        imps(1,:,:) = chebfun.jumpVals(funs);

        % Update impulsess:
        if ( size(imps, 1) > 1 )
           imps = [imps(1,:,:) ; newImps ; imps(2:end,:,:)];
        elseif any(newImps)
           imps(2,:,:) = newImps;
        end

    end
end

% Reassign data to f:
f.funs = funs;
f.funs = funs;
f.impulses = imps;

end