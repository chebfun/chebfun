function F = mtimes(F, G)
%*   MTIMES for SPHEREFUNV.
%   c*F or F*c multiplies each component of a SPHEREFUNV by a scalar.
%
%   A*F multiplies the vector of functions F by the matrix A assuming that
%   size(A,2) == size(F,1).
%
%   F*G calculates the inner product between F and G if size(F,3) ==
%   size(G,1). If the sizes are appropriate then F*G = dot(F.',G).
%
% See also SPHEREFUNV/TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) || isempty(G) )
    F = spherefunv;
    return
end

% If the SPHEREFUNV object is transposed, then compute (G.'*f.').'
if ( isa(F, 'spherefunv') && ~isa(G,  'spherefunv') )
    if ( F.isTransposed )
        F = mtimes(G.', F.');
        return
    end
end

if ( isa(G, 'spherefunv') && ~isa(F, 'spherefunv') )
    if ( G.isTransposed )
        F = mtimes( G.' , F.' ).' ;
        return
    end
end

if ( isa(F, 'double') )      % doubles * SPHEREFUNV
    if ( numel(F) == 1 )     % scalar * SPHEREFUNV
        const = F;
        F = G;
        for j = 1:3
            F.components{j} = const * F.components{j};
        end
    elseif ( size(F, 2) == 3 )   % matrix * column SPHEREFUNV
        vec = F;
        if ( size(vec, 1) == 1 ) 
            F = vec(1, 1) * G.components{1};
            for jj = 2:3
                F = F + vec(1, jj) * G.components{jj};
            end
        else
            store = {}; 
            for jj = 1:size(vec, 1) 
                store{jj} = mtimes(vec(jj,:), G); 
            end
            F = spherefun(store); 
        end
    else
        error('SPHEREFUN:SPHEREFUNV:mtimes:double1', 'Dimension mismatch.');
    end
    
elseif( isa(G, 'double') )        % SPHEREFUNV * double
    
    if ( numel(G) == 1 )          % SPHEREFUNV * scalar
        F = mtimes(G, F);
    else
        error('SPHEREFUN:SPHEREFUNV:mtimes:double2', ...
            'SPHEREFUNV and double size mismatch.');
    end
elseif ( isa(F,'spherefunv') && isa(G,'spherefunv') ) 
    
    % Dot product if dimensions are right:
    if ( F.isTransposed && ~G.isTransposed )
        F = dot(F, G);
    else
        error('SPHEREFUN:SPHEREFUNV:mtimes:sizes', 'Dimensions mismatch.');
    end
    
elseif ( isa(F,'spherefunv') && isa(G,'spherefun') )
    
    F = mtimes(G , F);
    
else 
    error('SPHEREFUN:SPHEREFUNV:mtimes:inputs', ...
        'SPHEREFUNV can only mtimes to SPHEREFUNV or double');
end

end
