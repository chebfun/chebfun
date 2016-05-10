function f = mtimes(f, g)
%*   CHEBFUN multiplication.
%   A*F and F*A multiplies the CHEBFUN F by the scalar A.
%
%   If F is an m-by-Inf row CHEBFUN and G is an Inf-by-n column CHEBFUN, F*G
%   returns the m-by-n matrix of pairwise inner products. F and G must have
%   the same domain.
%
% See also TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsChebfun = isa(f, 'chebfun');
gIsChebfun = isa(g, 'chebfun');

if ( isempty(f) || isempty(g) )             % [] * CHEBFUN or CHEBFUN * []

    f = [];

elseif ( fIsChebfun && gIsChebfun )         % CHEBFUN * CHEBFUN

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
        % Compute the inner product (we call CONJ() here because INNERPRODUCT()
        % is semilinear in the first factor, and we want to undo that):
        f = innerProduct(conj(f), g);
    else                                           % Column times row.
        % Outer-product of two CHEBFUNs is a CHEBFUN2.
        f = chebfun2.outerProduct(f, g);
    end

elseif ( fIsChebfun && isnumeric(g) )       % CHEBFUN * double

    if ( isscalar(g) )
        f = times(f, g);
    else

        % Check the dimensions.
        if ( (size(g, 1) ~= size(f, 2)) || (ndims(g) > 2) )
            error('CHEBFUN:CHEBFUN:mtimes:dims', ...
                'Matrix dimensions must agree.');
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
    end

elseif ( isnumeric(f) && gIsChebfun )       % double * CHEBFUN

        f = mtimes(g.', f.').';

elseif ( fIsChebfun )                       % CHEBFUN * ??? 

        f = (g.'*f.').';

else                                        % ??? * CHEBFUN or CHEBFUN * ???

    error('CHEBFUN:CHEBFUN:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));

end

end
