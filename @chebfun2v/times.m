function F = times( F , G ) 
%.* Times of two CHEBFUN2V objects. 
%   F.*G if F is a CHEBFUN2V and G is double returns the CHEBFUN2V after
%   componentwise multiplication.
%
%   F.*G if F is a double and G is a CHEBFUN2V returns the CHEBFUN2V after
%   componentwise multiplication.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    F = chebfun2v;
    return
end

if ( ~isa(F, 'chebfun2v') ) 
    F = times(G, F); 
    return
end

nF = F.nComponents; 

if ( isa(G, 'double') )             % CHEBFUN2V.*double
    if ( numel(G) == 1 )            % CHEBFUN2V.*scalar
        scalar = G;
        for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, scalar); 
        end
    elseif ( (size(G, 1) == nF) || ( F.isTransposed && (size(G, 2) == nF) ) )
        for jj = 1 : nF 
            F.components{jj} = times( F.components{jj}, G(jj) ); 
        end   
    else
        error('CHEBFUN:CHEBFUN2V:times:double', ...
            'CHEBFUN2V and double size mismatch.');
    end  
    
elseif ( isa(G, 'chebfun2v') )      % CHEBFUN2V . * CHEBFUN2V
    nG = G.nComponents; 
    if ( nF ~= nG ) 
         error('CHEBFUN:CHEBFUN2V:times:times', ...
             'CHEBFUN2V components mismatch.');
    end
    for jj = 1:nF 
        F.components{jj} = times(F.components{jj}, G.components{jj}); 
    end
    
elseif ( isa(G, 'chebfun2') )       % CHEBFUN2 * CHEBFUN2V
    for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, G); 
    end
    
else  % error
    error( 'CHEBFUN:CHEBFUN2V:times:inputs', 'Unrecognized input arguments.' );
end

end

