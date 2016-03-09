function F = times( F , G ) 
%.* Times of two DISKFUNV objects. 
%   F.*G if F is a DISKFUNV and G is double returns the DISKFUNV after
%   componentwise multiplication.
%
%   F.*G if F is a double and G is a DISKFUNV returns the DISKFUNV after
%   componentwise multiplication.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    F = diskfunv;
    return
end

if ( ~isa(F, 'diskfunv') ) 
    F = times(G, F); 
    return
end

nF = F.nComponents; 

if ( isa(G, 'double') )             % DISKFUNV.*double
    if ( numel(G) == 1 )            % DISKFUNV.*scalar
        scalar = G;
        for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, scalar); 
        end
    elseif ( (size(G, 1) == nF) || ( F.isTransposed && (size(G, 2) == nF) ) )
        for jj = 1 : nF 
            F.components{jj} = times( F.components{jj}, G(jj) ); 
        end   
    else
        error('CHEBFUN:DISKFUNV:times:double', ...
            'DISKFUNV and double size mismatch.');
    end  
    
elseif ( isa(G, 'diskfunv') )      % DISKFUNV . * DISKFUNV
    nG = G.nComponents; 
    if ( nF ~= nG ) 
         error('CHEBFUN:DISKFUNV:times:times', ...
             'DISKFUNV components mismatch.');
    end
    for jj = 1:nF 
        F.components{jj} = times(F.components{jj}, G.components{jj}); 
    end
    
elseif ( isa(G, 'diskfun') )       % DISKFUN * DISKFUNV
    for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, G); 
    end
    
else  % error
    error( 'CHEBFUN:DISKFUNV:times:inputs', 'Unrecognized input arguments.' );
end

end

