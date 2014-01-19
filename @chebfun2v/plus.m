function F = plus( F, G ) 
% + PLUS of two chebfun2v objects. 
% 
% F + G if F and G are chebfun2v objects does componentwise addition. 
% 
% F + G if F is a double and G is a chebfun2v does componentwise addition. 
% 
% F + G if F is a chebfun2v and G is a double does componentwise addition. 
% 
% plus(F,G) is called for the syntax F + G. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) ) 
    F = chebfun2v; 
    return
end

if ( ~isa( F , 'chebfun2v' ) )
    F = plus(G, F); 
    return
end

nF = F.nComponents;

if ( isa(G, 'double') ) 
    if ( numel(G) == 1 )
        for jj = 1 : nF 
            F.components{jj} = plus(F.components{jj}, G);
        end
    elseif ( numel(G) == nF ) 
        for jj = 1 : nF 
             F.components{jj} = plus(F.components{jj}, G(jj));
        end          
    else
        error('CHEBFUN2V:PLUS:doubleSize', 'Dimension mismatch.')
    end
elseif ( isa(G, 'chebfun2') ) 
    for jj = 1 : nF 
        F.components{jj} = plus(F.components{jj}, G);
    end
elseif ( isa(G, 'chebfun2v') )
    nG = G.nComponents; 
    if ( nG ~= nF ) 
        error('CHEBFUN2V:PLUS:components', 'The chebfun2v objects do not have the same components.')
    end
    if ( G.isTransposed ~= F.isTransposed )
        error('CHEBFUN2V:PLUS:Transposed', 'Dimension mismatch.')
    end
    for jj = 1 : nF 
        F.components{jj} = plus(F.components{jj}, G.components{jj});
    end
else
    error('CHEBFUN2V:PLUS:Type', 'Unrecongized input arguments')
end

end
    