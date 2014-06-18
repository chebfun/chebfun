function f = mtimes(f, g)
%*   CHEBFUN multiplication.
%   A*F and F*A multiplies the CHEBFUN F by the scalar A.
%
%   If F is an m-by-Inf row CHEBFUN and G is an Inf-by-n column CHEBFUN, F*G
%   returns the m-by-n matrix of pairwise inner products. F and G must have
%   the same domain.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(f, 'chebfun') )   % ??? * CHEBFUN

    % Ensure CHEBFUN is the first input:
    if ( ~g(1).isTransposed )
        f = mtimes(g, f);
    else
        f = mtimes(g.', f.').';
    end

elseif ( isempty(g) )       % CHEBFUN * []

    f = [];
    
elseif ( isnumeric(g) )     % CHEBFUN * double

    if ( isscalar(g) )
        f = times(f, g);
        return
    end

    if ( numel(f) == 1 )
        % Array-valued CHEBFUN case:
    
        % Loop over the FUNs:
        for k = 1:numel(f.funs)
            f.funs{k} = mtimes(f.funs{k}, g);
        end

        % Multiply the pointValues:
        f.pointValues = f.pointValues * g;

    else
        % QUASIMATRIX case:
        numCols = numel(f);
        if ( length(g) ~= numCols && min(size(g)) ~= 1 )
            error('CHEBFUN:CHEBFUN:mtimes:dims', ...
                'Matrix dimensions must agree.');
        end
        % Transpose g if f is a row CHEBFUN:
        if ( f(1).isTransposed )
            g = g.';
        end
        s = f(1)*g(1,:);
        % Loop over the columns:
        for k = 2:numCols
            s = s + f(k).*g(k,:);
        end
        f = s;
    end

elseif ( ~isa(g, 'chebfun') )

    error('CHEBFUN:CHEBFUN:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
else                        % CHEBFUN' * CHEBFUN

    % We can't do MTIMES() on two CHEBFUNs that have the same transpose state.
    if ( f(1).isTransposed == g(1).isTransposed )
        if ( numColumns(f) ~= numColumns(g) )
            error('CHEBFUN:CHEBFUN:mtimes:dims', ...
                ['Matrix dimensions must agree. Use f.*g to multiply ' ...
                 'two CHEBFUN objects.']);
        else
            error('CHEBFUN:CHEBFUN:mtimes:dims', ...
                'Matrix dimensions must agree.');
        end
    end

    if ( f(1).isTransposed && ~g(1).isTransposed ) % Row times column.
        
        f = innerProduct(conj(f), g);
%         % Overlap:
%         [f, g] = overlap(f, g);
% 
%         % Compute the inner product (we call CONJ() here because INNERPRODUCT()
%         % is semilinear in the first factor, and we want to undo that):
%         S = 0;
%         for k = 1:numel(f.funs)
%             S = S + innerProduct(conj(f.funs{k}), g.funs{k});
%         end
% 
%         % Output in f:
%         f = S;
    else                                     % Column times row.
        % Outer-product of two chebfuns is a chebfun2. 
        f = chebfun2.outerProduct(f, g);
    end

end

end
