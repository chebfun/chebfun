function F = plus( F, G ) 
%+   PLUS of two SPHEREFUNV objects. 
%   F + G if F and G are SPHEREFUNV objects does componentwise addition. 
%
%   F + G if F is a double and G is a SPHEREFUNV does componentwise
%   addition.
% 
%   F + G if F is a SPHEREFUNV and G is a double does componentwise
%   addition.
% 
%   PLUS(F, G) is called for the syntax 'F + G'. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty(F) || isempty(G) ) 
    F = spherefunv; 
    return
end

if ( ~isa(F , 'spherefunv') )
    F = plus(G, F); 
    return
end

if ( isa(G, 'double') )              % SPHEREFUNV + DOUBLE
    if ( numel(G) == 1 )             % SPHEREFUNV + SCALAR
        for jj = 1 : 3 
            F.components{jj} = plus(F.components{jj}, G);
        end
    elseif ( numel(G) == 3 )        % SPHEREFUNV + MATRIX
        for jj = 1 : 3 
             F.components{jj} = plus(F.components{jj}, G(jj));
        end          
    else
        error('SPHEREFUN:SPHEREFUNV:plus:doubleSize', 'Dimension mismatch.')
    end
elseif ( isa(G, 'spherefun') )        % SPHEREFUNV + SPHEREFUN
    for jj = 1:3 
        F.components{jj} = plus(F.components{jj}, G);
    end
elseif ( isa(G, 'spherefunv') )       % SPHEREFUNV + SPHEREFUNV
    if ( G.isTransposed ~= F.isTransposed )
        error('SPHEREFUN:SPHEREFUNV:plus:transposed', 'Dimension mismatch.')
    end
    for jj = 1:3                  % Add each component together
        F.components{jj} = plus(F.components{jj}, G.components{jj});
    end
else
    error('SPHEREFUN:SPHEREFUNV:plus:type', 'Unrecongized input arguments.')
end

end
    
