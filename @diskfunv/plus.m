function F = plus( F, G ) 
% + PLUS of two DISKFUNV objects. 
%   F + G if F and G are DISKFUNV objects does componentwise addition. 
%
%   F + G if F is a double and G is a DISKFUNV does componentwise addition. 
% 
%   F + G if F is a DISKFUNV and G is a double does componentwise addition.
% 
%   PLUS(F,G) is called for the syntax F + G. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) ) 
    F = diskfunv; 
    return
end

if ( ~isa( F , 'diskfunv' ) )
    F = plus(G, F); 
    return
end

nF = F.nComponents;

if ( isa(G, 'double') )              % DISKFUNV + DOUBLE
    if ( numel(G) == 1 )             % DISKFUNV + SCALAR
        for jj = 1 : nF 
            F.components{jj} = plus(F.components{jj}, G);
        end
    elseif ( numel(G) == nF )        % DISKFUNV + MATRIX
        for jj = 1 : nF 
             F.components{jj} = plus(F.components{jj}, G(jj));
        end          
    else
        error('CHEBFUN:DISKFUNV:plus:doubleSize', 'Dimension mismatch.')
    end
elseif ( isa(G, 'diskfun') )        % DISKFUNV + CHEBFUN
    for jj = 1 : nF 
        F.components{jj} = plus(F.components{jj}, G);
    end
elseif ( isa(G, 'diskfunv') )       % DISKFUNV + DISKFUNV
    nG = G.nComponents; 
    if ( nG ~= nF ) 
        error('CHEBFUN:DISKFUNV:plus:components', ...
            'The DISKFUNV objects do not have the same components.')
    end
    if ( G.isTransposed ~= F.isTransposed )
        error('CHEBFUN:DISKFUNV:plus:transposed', 'Dimension mismatch.')
    end
    for jj = 1 : nF                  % Add each component together
        F.components{jj} = plus(F.components{jj}, G.components{jj});
    end
else
    error('CHEBFUN:DISKFUNV:plus:type', 'Unrecongized input arguments')
end

end
    
