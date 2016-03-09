function F = power( F, G )
%.^ Componentwise power for DISKFUNV.
%   F.^G where F is a DISKFUNV and G is a double returns the result from
%   componentwise powers.
%
%   F.^G where F is a double and G is a CHEBFUN2 returns from componentwise
%   powers.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) || isempty( G ) )
    F = diskfunv;
    return
end

if ( isa( F ,'double' ) )       % scalar . ^ DISKFUNV
    
    if ( numel(F) == 1 )
        scalar = F;
        F = G;
        for jj = 1 : F.nComponents
            F.components{jj} = power( scalar, G.components{jj} );
        end
    else
        error('CHEBFUN:DISKFUNV:power:double', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )      % DISKFUNV . ^ scalar
    
    if ( numel(G) == 1 )
        scalar = G;
        for jj = 1 : F.nComponents
            F.components{jj} = power( F.components{jj}, scalar );
        end
    else
        error('CHEBFUN:DISKFUNV:power:double', ...
            'DISKFUNV and double size mismatch.');
    end
    
elseif (isa(F,'diskfunv') && isa(G,'diskfunv') )  % DISKFUNV.^DISKFUNV
    
    error('CHEBFUN:DISKFUNV:power:size', ...
        'DISKFUNV dimension mismatch.');
    
else  % error
    
    error('CHEBFUN:DISKFUNV:power:inputs', ...
        'Unrecognized input arguments.');
    
end

end
