function F = rdivide( F, G )
%./   Pointwise DISKFUNV right divide.
%
% See also LDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
   F = diskfunv;
   return 
end

% Componentwise divide. 
if ( isa(F, 'diskfunv') && isa(G, 'diskfunv') )
    
    nF = F.nComponents; 
    nG = F.nComponents; 
    if ( nF ~= nG ) 
        error('CHEBFUN:DISKFUNV:rdivide:dim','Dimension mismatch.');
    end
    for jj = 1 : nF 
        F.components{jj} = rdivide( F.components{jj}, G.components{jj} ); 
    end
    
elseif ( isa(F, 'diskfunv') && ( isa(G, 'diskfun') || isa(G, 'double') ) ) 
    
    nF = F.nComponents; 
    for jj = 1 : nF 
        F.components{jj} = rdivide( F.components{jj}, G ); 
    end
    
elseif  ( isa(G, 'diskfunv') && ( isa(F, 'diskfun') || isa(F, 'double') ) ) 
    
    H = G; 
    nG = G.nComponents;
    for jj = 1 : nG
        H.components{jj} = rdivide( F, H.components{jj} ); 
    end
    F = H; 
    
else
    
   error('CHEBFUN:DISKFUNV:rdivide:inputs','Unrecognized input arguments.')
   
end
    
end
