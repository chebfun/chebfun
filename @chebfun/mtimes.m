function f = mtimes(f, g)
%*   CHEBFUN multiplication.
%   F.*G multiplies the CHEBFUN objects F and G or a CHEBFUN by a scalar if
%   either F or G is a scalar.

if ( ~isa(f, 'chebfun') )   % ??? * CHEBFUN

    % Ensure CHEBFUN is the first input:
    f = times(g, f);

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

    error('CHEBFUN:times:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type '] ...
          ['%s and %s.'], class(f), class(g));
else                        % CHEBFUN' * CHEBFUN

    if ( size(f, 1) ~= size(g, 2) )
        if ( (~isTransposed(f) && (size(f, 2) == size(g, 2))) || ...
             (isTransposed(f) && (size(f, 1) == size(g, 1))) )
            error('CHEBFUN:times:dims', ...
                ['Matrix dimensions must agree. Use f.*g to multiply ' ...
                 'two chebfun objects.']);
        else
            error('CHEBFUN:times:dims', ...
                'Matrix dimensions must agree.');
        end
    end

    % Overlap:
    [f, g] = overlap(f, g);
    
    % Compute the inner product:
    S = 0;
    for k = 1:numel(f.funs)
        S = S + innerProduct(f.funs{k}, g.funs{k});
    end
    
    % Output in f:
    f = S;

end

end
