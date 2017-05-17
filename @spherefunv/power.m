function F = power(F, G)
%.^   Componentwise power for SPHEREFUNV.
%   F.^G where F is a SPHEREFUNV and G is a double returns the result from
%   componentwise powers.
%
%   F.^G where F is a double and G is a SPHEREFUN returns from componentwise
%   powers.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) || isempty(G) )
    F = spherefunv;
    return
end

if ( isa(F ,'double') )       % scalar . ^ SPHEREFUNV
    
    if ( numel(F) == 1 )
        scalar = F;
        F = G;
        for jj = 1 : 3
            F.components{jj} = power( scalar, G.components{jj} );
        end
    else
        error('SPHEREFUN:SPHEREFUNV:power:double', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )    % SPHEREFUNV . ^ scalar
    
    if ( numel(G) == 1 )
        scalar = G;
        for jj = 1:3
            F.components{jj} = power(F.components{jj}, scalar);
        end
    else
        error('SPHEREFUN:SPHEREFUNV:power:double', ...
            'SPHEREFUNV and double size mismatch.');
    end
    
elseif (isa(F,'spherefunv') && isa(G,'spherefunv') )  % SPHEREFUNV.^SPHEREFUNV
    
    error('SPHEREFUN:SPHEREFUNV:power:size', ...
        'SPHEREFUNV dimension mismatch.');
    
else  % error
    
    error('SPHEREFUN:SPHEREFUNV:power:inputs', ...
        'Unrecognized input arguments.');
    
end

end
