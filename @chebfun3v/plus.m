function F = plus(F, G) 
% + PLUS of two CHEBFUN3V objects. 
%
%   F + G if F and G are CHEBFUN3V objects does componentwise addition. 
%
%   F + G if F is a double and G is a CHEBFUN3V does componentwise addition. 
% 
%   F + G if F is a CHEBFUN3V and G is a double does componentwise addition.
% 
%   PLUS(F,G) is called for the syntax F + G. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) || isempty( G ) ) 
    F = chebfun3v; 
    return
end

if ( ~isa( F , 'chebfun3v' ) )
    F = plus(G, F); 
    return
end

nF = F.nComponents;

if ( isa(G, 'double') )              % CHEBFUN3V + DOUBLE
    if ( numel(G) == 1 )             % CHEBFUN3V + SCALAR
        for jj = 1 : nF 
            F.components{jj} = plus(F.components{jj}, G);
        end
    elseif ( numel(G) == nF )        % CHEBFUN3V + VECTOR
        for jj = 1 : nF 
             F.components{jj} = plus(F.components{jj}, G(jj));
        end          
    else
        error('CHEBFUN:CHEBFUN3V:plus:doubleSize', 'Dimension mismatch.')
    end
elseif ( isa(G, 'chebfun3') )        % CHEBFUN3V + CHEBFUN3
    for jj = 1 : nF 
        F.components{jj} = plus(F.components{jj}, G);
    end
elseif ( isa(G, 'chebfun3v') )       % CHEBFUN3V + CHEBFUN3V
    nG = G.nComponents; 
    if ( nG ~= nF ) 
        error('CHEBFUN:CHEBFUN3V:plus:components', ...
            'The chebfun3v objects do not have the same components.')
    end
    if ( G.isTransposed ~= F.isTransposed )
        error('CHEBFUN:CHEBFUN3V:plus:transposed', 'Dimension mismatch.')
    end
    for jj = 1 : nF                  % Add each component together
        F.components{jj} = plus(F.components{jj}, G.components{jj});
    end
else
    error('CHEBFUN:CHEBFUN3V:plus:type', 'Unrecongized input arguments')
end

end
    
