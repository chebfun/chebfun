function F = power( F, G )
%.^ Componentwise power for chebfun2v.
%   F.^G where F is a chebfun2v and G is a double returns the result from
%   componentwise powers.
%
%   F.^G where F is a double and G is a chebfun2 returns from componentwise
%   powers.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check:
if ( isempty( F ) || isempty( G ) )
    F = chebfun2v;
    return
end

if ( isa( F ,'double' ) )       % scalar . ^ chebfun2v
    
    if ( numel(F) == 1 )
        scalar = F;
        F = G;
        for jj = 1 : F.nComponents
            F.components{jj} = power( scalar, G.components{jj} );
        end
    else
        error('CHEBFUN2v:mtimes:double', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )      % chebfun2v . ^ scalar
    
    if ( numel(G) == 1 )
        scalar = G;
        for jj = 1 : F.nComponents
            F.components{jj} = power( F.components{jj}, scalar );
        end
    else
        error('CHEBFUN2v:mtimes:double', 'Chebfun2v and double size mismatch.');
    end
    
elseif (isa(F,'chebfun2v') && isa(G,'chebfun2v') )  % chebfun2v.^chebfun2v
    
    error('CHEBFUN2v:power:size', 'Chebfun2v dimension mismatch.');
    
else  % error
    
    error('CHEBFUN2v:power:inputs', 'Unrecognized input arguments.');
    
end

end