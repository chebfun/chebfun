function F = mtimes( F, G )
%*  mtimes for DISKFUNV.
%
%  c*F or F*c multiplies each component of the DISKFUNV F by the scalar c.
%
%  A*F multiplies the DISKFUNV F by the matrix A assuming that
%  size(A,2) == size(F,1).
%
%  F*G calculates the inner product between F and G if size(F,3) ==
%  size(G,1). If the sizes are compatible then F*G = dot(F.',G).
%
% See also TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( ( isempty(F) ) || ( isempty(G) ) )
    F = diskfunv;
    return
end

% If the DISKFUNV object is transposed, then compute (G.'*f.').':
if ( isa( F, 'diskfunv' ) && ~isa( G,  'diskfunv' ) )
    if ( F.isTransposed )
        F = mtimes( G.', F.' );
        return
    end
end

% Transpose if G is the only DISKFUNV input:  
if ( isa(G, 'diskfunv') && ~isa(F, 'diskfunv') )
    if ( G.isTransposed )
        F = mtimes( G.' , F.' ).' ;
        return
    end
end

if ( isa( F, 'double' ) )      % doubles * DISKFUNV
    if ( numel(F) == 1 )       % scalar * DISKFUNV
        const = F;
        F = G;
        for j = 1:F.nComponents
            F.components{j} = const * F.components{j};
        end
    elseif ( size(F, 2) == G.nComponents )   % matrix * column DISKFUNV
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
            F = diskfunv(store); 
        end
    else
        error('DISKFUN:diskfunv:mtimes:double1', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )          % DISKFUNV * double
    
    if ( numel( G ) == 1 )          % DISKFUNV * scalar
        F = mtimes( G, F );
    else
        error('DISKFUN:DISKFUNV:mtimes:double2', ...
            'DISKFUNV and double size mismatch.');
    end
elseif ( isa(F, 'diskfunv') && isa(G, 'diskfunv') ) % dot product if dimensions are riGht.
    
    if ( ( F.isTransposed ) && ( ~G.isTransposed ) )
        F = dot( F, G );
    else
        error('DISKFUN:diskfunv:mtimes:sizes', 'Dimensions mismatch.');
    end
    
elseif ( isa(F, 'diskfunv') && isa(G, 'diskfun') )
    
    F = mtimes( G , F );
    
else 
    error('DISKFUN:DISKFUNV:mtimes:inputs', ...
        'DISKFUNV can only mtimes to DISKFUNV or double');
end
end