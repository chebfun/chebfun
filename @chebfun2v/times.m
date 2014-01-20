function F = times( F , G ) 
%.* Times of two chebfun2v objects. 
%
% F.*G if F is a chebfun2v and G is double returns the chebfun2v after
% componentwise multiplication. 
%
% F.*G if F is a double and G is a chebfun2v returns the chebfun2v after
% componentwise multiplication.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

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

if ( isa(G, 'double') )            % chebfun2v.*double
    if ( numel(G) == 1 )           % chebfun2v.*scalar
        scalar = G;
        for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, scalar); 
        end
    elseif ( ( size(G, 1) == nF ) || ( ( F.isTransposed ) && ( size(G, 2) == nF )) )
        for jj = 1 : nF 
            F.components{jj} = times( F.components{jj}, G(jj) ); 
        end   
    else
            error('CHEBFUN2v:times:double','Chebfun2v and double size mismatch.');
    end  
elseif ( isa(G, 'chebfun2v') )      % chebfun2v . * chebfun2v
    nG = G.nComponents; 
    if ( nF ~= nG ) 
         error('CHEBFUN2v:times','Chebfun2v components mismatch.');
    end
    for jj = 1:nF 
        F.components{jj} = times(F.components{jj}, G.components{jj}); 
    end
    
elseif ( isa(G, 'chebfun2') )        % chebfun2 * chebfun2v
    
    for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, G); 
    end
    
else  % error
    error( 'CHEBFUN2v:times:inputs', 'Unrecognized input arguments.' );
end
end