function F = rdivide( F, G )
%./   Pointwise CHEBFUN2V right divide.
%
% See also LDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
   F = chebfun2v;
   return 
end

% Componentwise divide. 
if ( isa(F, 'chebfun2v') && isa(G, 'chebfun2v') )
    
    nF = F.nComponents; 
    nG = F.nComponents; 
    if ( nF ~= nG ) 
        error('CHEBFUN:CHEBFUN2V:rdivide:dim','Dimension mismatch.');
    end
    for jj = 1 : nF 
        F.components{jj} = rdivide( F.components{jj}, G.components{jj} ); 
    end
    
elseif ( isa(F, 'chebfun2v') && ( isa(G, 'chebfun2') || isa(G, 'double') ) ) 
    
    nF = F.nComponents; 
    for jj = 1 : nF 
        F.components{jj} = rdivide( F.components{jj}, G ); 
    end
    
elseif  ( isa(G, 'chebfun2v') && ( isa(F, 'chebfun2') || isa(F, 'double') ) ) 
    
    H = G; 
    nG = G.nComponents;
    for jj = 1 : nG
        H.components{jj} = rdivide( F, H.components{jj} ); 
    end
    F = H; 
    
else
    
   error('CHEBFUN:CHEBFUN2V:rdivide:inputs','Unrecognized input arguments.')
   
end
    
end
