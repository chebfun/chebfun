function F = power(F, G)
%.^   Componentwise power for CHEBFUN3V objects.
%   F.^G computes componentwise power of a CHEBFUN3V object F and a double
%   G.
%
%   Alternatively, F can be a double if G is a CHEBFUN3V object.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) || isempty(G) )
    F = chebfun3v();
    return
end

if ( isa(F ,'double') )                              % scalar .^ CHEBFUN3V
    
    if ( numel(F) == 1 )
        scalar = F;
        F = G;
        for jj = 1:F.nComponents
            F.components{jj} = power(scalar, G.components{jj});
        end
    else
        error('CHEBFUN:CHEBFUN3V:power:double', 'Dimension mismatch.');
    end
    
elseif ( isa(G, 'double') )                          % CHEBFUN3V .^ scalar
    
    if ( numel(G) == 1 )
        scalar = G;
        for jj = 1:F.nComponents
            F.components{jj} = power(F.components{jj}, scalar);
        end
    else
        error('CHEBFUN:CHEBFUN3V:power:double', ...
            'CHEBFUN3V and double size mismatch.');
    end
    
elseif (isa(F,'chebfun3v') && isa(G,'chebfun3v') )  % CHEBFUN3V.^CHEBFUN3V
    error('CHEBFUN:CHEBFUN3V:power:size', ...
        'CHEBFUN3V dimension mismatch.');
    
else  % error
    error('CHEBFUN:CHEBFUN3V:power:inputs', ...
        'Unrecognized input arguments.');
end

end