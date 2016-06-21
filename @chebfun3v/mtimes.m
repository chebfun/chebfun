function F = mtimes(F, G)
%*   MTIMES for CHEBFUN3V objects.
%   c*F or F*c multiplies each component of a CHEBFUN3V object F by a 
%   scalar c.
%
%   A*F multiplies the vector of CHEBFUN3 objects F by the matrix A 
%   assuming that size(A, 2) == size(F, 1). Note that size(A, 1) must be 
%   at most 3.
%
%   F*G calculates the inner product between F and G if size(F, 4) ==
%   size(G, 1). If the sizes are appropriate then F*G = dot(F.', G).
%
% See also CHEBFUN3V/TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( ( isempty(F) ) || ( isempty(G) ) )
    F = chebfun3v;
    return
end

% If the input CHEBFUN3V object is transposed, then make the output also a
% row CHEBFUN3V object:
if ( isa(F, 'chebfun3v') && ~isa(G, 'chebfun3v') )
    if ( F.isTransposed )
        F = mtimes(G.', F.').';
        return
    end
end

if ( isa(G, 'chebfun3v') && ~isa(F, 'chebfun3v') )
    if ( G.isTransposed )
        F = mtimes(G.', F.').';
        return
    end
end

if ( isa(F, 'double') )        % doubles * CHEBFUN3V
    if ( numel(F) == 1 )       % scalar * CHEBFUN3V
        const = F;
        F = G;
        for j = 1:F.nComponents
            F.components{j} = const * F.components{j};
        end
    elseif ( size(F, 2) == G.nComponents )   % matrix * column CHEBFUN3V
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
            F = chebfun3v(store);
        end
    else
        error('CHEBFUN:CHEBFUN3V:mtimes:double', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )          % CHEBFUN3V * double
    
    if ( numel(G) == 1 )          % CHEBFUN3V * scalar
        F = mtimes(G, F);
    else
        error('CHEBFUN:CHEBFUN3V:mtimes:double', ...
            'CHEBFUN3V and double size mismatch.');
    end
    
elseif ( isa(F, 'chebfun3v') && isa(G, 'chebfun3v') ) 
    % dot product if dimensions are right.
    
    if ( F.isTransposed && ~G.isTransposed )
        F = dot(F, G);
    else
        error('CHEBFUN:CHEBFUN3V:mtimes:sizes', 'Dimensions mismatch.');
    end
    
elseif ( isa(F,'chebfun3v') && isa(G,'chebfun3') )
    F = mtimes(G , F);
    
else
    error('CHEBFUN:CHEBFUN3V:mtimes:inputs', ...
        'CHEBFUN3V can only mtimes to CHEBFUN3V or double');
end

end