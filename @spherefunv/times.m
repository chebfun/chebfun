function F = times( F , G ) 
%.*   Times of two SPHEREFUNV objects. 
%   F.*G if F is a SPHEREFUNV and G is double returns the SPHEREFUNV after
%   componentwise multiplication.
%
%   F.*G if F is a SPHEREFUN and G is a SPHEREFUNV returns the SPHEREFUNV
%   after multiplication of F by each component of G.
%
%   F.*G if F is a double and G is a SPHEREFUNV returns the SPHEREFUNV
%   after multiplication of F by each component of G.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
    F = spherefunv;
    return
end

if ( ~isa(F, 'spherefunv') ) 
    F = times(G, F); 
    return
end

nF = 3; 

if ( isa(G, 'double') )             % SPHEREFUNV.*double
    if ( numel(G) == 1 )            % SPHEREFUNV.*scalar
        scalar = G;
        for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, scalar); 
        end
    elseif ( ( size(G, 1) == nF ) || ( F.isTransposed && (size(G, 2) == nF) ) )
        for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, G(jj)); 
        end   
    else
        error('SPHEREFUN:SPHEREFUNV:times:double', ...
            'SPHEREFUNV and double size mismatch.');
    end  
    
elseif ( isa(G, 'spherefunv') )      % SPHEREFUNV . * SPHEREFUNV
    for jj = 1:3 
        F.components{jj} = times(F.components{jj}, G.components{jj}); 
    end
    
elseif ( isa(G, 'spherefun') )       % SPHEREFUN * SPHEREFUNV
    for jj = 1 : nF 
        F.components{jj} = times(F.components{jj}, G); 
    end
    
else  % error
    error( 'SPHEREFUN:SPHEREFUNV:times:inputs', ...
        'Unrecognized input arguments.');
end

end
