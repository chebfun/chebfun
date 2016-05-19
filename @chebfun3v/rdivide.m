function F = rdivide(F, G)
%./   Pointwise right divide for CHEBFUN3V objects.
%
% See also CHEBFUN3V/LDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
   F = chebfun3v;
   return 
end

% Componentwise divide. 
if ( isa(F, 'chebfun3v') && isa(G, 'chebfun3v') )
    nF = F.nComponents; 
    nG = F.nComponents; 
    if ( nF ~= nG ) 
        error('CHEBFUN:CHEBFUN3V:rdivide:dim','Dimension mismatch.');
    end
    for jj = 1:nF
        F.components{jj} = rdivide(F.components{jj}, G.components{jj});
    end
    
elseif ( isa(F, 'chebfun3v') && ( isa(G, 'chebfun3') || isa(G, 'double') ) ) 
    nF = F.nComponents; 
    for jj = 1:nF
        F.components{jj} = rdivide(F.components{jj}, G);
    end
    
elseif  ( isa(G, 'chebfun3v') && ( isa(F, 'chebfun3') || isa(F, 'double') ) )     
    H = G; 
    nG = G.nComponents;
    for jj = 1:nG
        H.components{jj} = rdivide(F, H.components{jj});
    end
    F = H;
    
else
   error('CHEBFUN:CHEBFUN3V:rdivide:inputs','Unrecognized input arguments.')
end
    
end