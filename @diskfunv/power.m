function F = power( F, G )
%.^ Componentwise power for DISKFUNV.
%   F.^G where F is a DISKFUNV and G is a double returns the result from
%   componentwise powers.
%
%   F.^G where F is a double and G is a CHEBFUN2 returns from componentwise
%   powers.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) || isempty( G ) )
    F = diskfunv;
    return
end

if ( isa( F, 'double' ) )       % scalar . ^ DISKFUNV
    
    if ( numel(F) == 1 )
        
        scalar = F;
        F = G;  
        F.components{1} = power( scalar, G.components{1} );
        F.components{2} = power( scalar, G.components{2} );
        
    else
        
        error('CHEBFUN:DISKFUNV:power:double', 'Dimension mismatch.');
    
    end
    
elseif ( isa(G, 'double') )      % DISKFUNV . ^ scalar
    
    if ( numel(G) == 1 )
        
        scalar = G;
        F.components{1} = power( F.components{1}, scalar );
        F.components{2} = power( F.components{2}, scalar );
        
    else
        
        error('CHEBFUN:DISKFUNV:power:double',...
                     'DISKFUNV and scalar dimension mismatch.');
                 
    end
    
elseif ( isa(F,'diskfunv') && isa(G,'diskfunv') )  % DISKFUNV.^DISKFUNV
    
    error('CHEBFUN:DISKFUNV:power:size', 'DISKFUNV dimension mismatch.');
    
else  % error
    
    error('CHEBFUN:DISKFUNV:power:inputs', 'Unrecognized input arguments.');
    
end

end
