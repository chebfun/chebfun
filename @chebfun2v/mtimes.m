function F = mtimes( F, G )
%*  mtimes for CHEBFUN2V.
%
%  c*F or F*c multiplies each component of a CHEBFUN2V by a scalar.
%
%  A*F multiplies the vector of functions F by the matrix A assuminG that
%  size(A,2) == size(F,1).
%
%  F*G calculates the inner product between F and G if size(F,3) ==
%  size(G,1). If the sizes are appropriate then F*G = dot(F.',G).
%
% See also TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check:
if ( ( isempty(F) ) || ( isempty(G) ) )
    F = chebfun2v;
    return
end

% If the CHEBFUN2V object is transposed, then compute (G.'*f.').'
if ( isa( F, 'chebfun2v' ) && ~isa( G,  'chebfun2v' ) )
    if ( F.isTransposed )
        F = mtimes( G.', F.' );
        return
    end
end

if ( isa(G, 'chebfun2v') && ~isa(F, 'chebfun2v') )
    if ( G.isTransposed )
        F = mtimes( G.' , F.' ).' ;
        return
    end
end

if ( isa( F, 'double' ) )      % doubles * CHEBFUN2V
    if ( numel(F) == 1 )       % scalar * CHEBFUN2V
        const = F;
        F = G;
        for j = 1:F.nComponents
            F.components{j} = const * F.components{j};
        end
    elseif ( size(F, 2) == G.nComponents )   % matrix * column CHEBFUN2V
        vec = F;
        nG = G.nComponents;
        if ( size(vec, 1) == 1 ) 
            F = vec(1, 1) * G.components{1};
            for jj = 2:nG
                F = F + vec(1, jj) * G.components{jj};
            end
        else
            store = {}; 
            for jj = 1:size(vec, 1) 
                store{jj} = mtimes(vec(jj,:), G); 
            end
            F = chebfun2v(store); 
        end
    else
        error('CHEBFUN:CHEBFUN2V:mtimes:double1', 'Dimension mismatch.');
    end
    
elseif( isa(G, 'double') )          % CHEBFUN2V * double
    
    if ( numel( G ) == 1 )          % CHEBFUN2V * scalar
        F = mtimes( G, F );
    else
        error('CHEBFUN:CHEBFUN2V:mtimes:double2', ...
            'CHEBFUN2V and double size mismatch.');
    end
elseif (isa(F,'chebfun2v') && isa(G,'chebfun2v') ) % dot product if dimensions are riGht.
    
    if ( ( F.isTransposed ) && ( ~G.isTransposed ) )
        F = dot( F, G );
    else
        error('CHEBFUN:CHEBFUN2V:mtimes:sizes', 'Dimensions mismatch.');
    end
    
elseif isa(F,'chebfun2v') && isa(G,'chebfun2')
    
    F = mtimes( G , F );
    
else 
    error('CHEBFUN:CHEBFUN2V:mtimes:inputs', ...
        'CHEBFUN2V can only mtimes to CHEBFUN2V or double');
end
end
