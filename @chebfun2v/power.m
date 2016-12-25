function F = power( F, G )
%.^ Componentwise power for CHEBFUN2V.
%   F.^G where F is a CHEBFUN2V and G is a double returns the result from
%   componentwise powers.
%
%   F.^G where F is a double and G is a CHEBFUN2 returns from componentwise
%   powers.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) || isempty( G ) )
    F = chebfun2v;
    return
end

if ( isa( F ,'double' ) )       % scalar . ^ CHEBFUN2V
    
    if ( numel(F) == 1 )
        scalar = F;
        F = G;
        for jj = 1 : F.nComponents
            F.components{jj} = power( scalar, G.components{jj} );
        end
    else
        error('CHEBFUN:CHEBFUN2V:power:double', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )      % CHEBFUN2V . ^ scalar
    
    if ( numel(G) == 1 )
        scalar = G;
        for jj = 1 : F.nComponents
            F.components{jj} = power( F.components{jj}, scalar );
        end
    else
        error('CHEBFUN:CHEBFUN2V:power:double', ...
            'CHEBFUN2V and double size mismatch.');
    end
    
elseif (isa(F,'chebfun2v') && isa(G,'chebfun2v') )  % CHEBFUN2V.^CHEBFUN2V
    
    error('CHEBFUN:CHEBFUN2V:power:size', ...
        'CHEBFUN2V dimension mismatch.');
    
else  % error
    
    error('CHEBFUN:CHEBFUN2V:power:inputs', ...
        'Unrecognized input arguments.');
    
end

end
