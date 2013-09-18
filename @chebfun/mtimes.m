function f = mtimes(f, g)
%*   CHEBFUN multiplication.
%   A*F and F*A multiplies the CHEBFUN F by the scalar A.
%
%   If F is an m-by-Inf row CHEBFUN and G is an Inf-by-n column CHEBFUN, F*G
%   returns the m-by-n matrix of pairwise inner products.  F and G must have
%   the same domain.
%
%   See also TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( ~isa(f, 'chebfun') )   % ??? * CHEBFUN

    % Ensure CHEBFUN is the first input:
    f = mtimes(g, f);

elseif ( isempty(g) )       % CHEBFUN * []

    f = [];
    
elseif ( isnumeric(g) )     % CHEBFUN * double

    % Loop over the FUNs:
    for k = 1:numel(f.funs)
        f.funs{k} = mtimes(f.funs{k}, g);
    end

    % Multiply the impulses:
    f.impulses = f.impulses * g;

elseif ( ~isa(g, 'chebfun') )

    error('CHEBFUN:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
else                        % CHEBFUN' * CHEBFUN

    % We can't do MTIMES() on two CHEBFUNs that have the same transpose state.
    if ( f.isTransposed == g.isTransposed )
        if ( ~all(size(f) == size(g)) )
            error('CHEBFUN:mtimes:dims', ...
                ['Matrix dimensions must agree. Use f.*g to multiply ' ...
                 'two chebfun objects.']);
        else
            error('CHEBFUN:mtimes:dims', ...
                'Matrix dimensions must agree.');
        end
    end

    if ( f.isTransposed && ~g.isTransposed ) % Row times column.
        % Overlap:
        [f, g] = overlap(f, g);

        % Compute the inner product (we call CONJ() here because INNERPRODUCT()
        % is semilinear in the first factor, and we want to undo that):
        S = 0;
        for k = 1:numel(f.funs)
            S = S + innerProduct(conj(f.funs{k}), g.funs{k});
        end

        % Output in f:
        f = S;
    else                                     % Column times row.
        % [TODO]:  Implement (and document) this once we have CHEBFUN2.
        error('CHEBFUN:mtimes:colTimesRow', ...
              'Support for (column)*(row) products not yet implemented.');
    end

end

end
