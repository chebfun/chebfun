function F = rdivide( F, G )
%./   Pointwise DISKFUNV right divide.
%
% See also LDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
   F = diskfunv();
   return 
end

% Componentwise divide: 
if ( isa(F, 'diskfunv') && isa(G, 'diskfunv') )  
    
    % DISKFUNV./DISKFUNV
    F.components{1} = rdivide( F.components{1}, G.components{1} ); 
    F.components{2} = rdivide( F.components{2}, G.components{2} ); 
     
elseif ( isa(F, 'diskfunv') && ( isa(G, 'diskfun') || isa(G, 'double') ) ) 
    
    % DISKFUNV./DISKFUN or DISKFUNV./DOUBLE:
    F.components{1} = rdivide( F.components{1}, G ); 
    F.components{2} = rdivide( F.components{2}, G ); 
       
elseif  ( isa(G, 'diskfunv') && ( isa(F, 'diskfun') || isa(F, 'double') ) ) 
    
    % DISKFUN./DISKFUNV or  DOUBLE./DISKFUNV:
    H = G;   
    H.components{1} = rdivide( F, H.components{1} ); 
    H.components{2} = rdivide( F, H.components{2} );
    F = H; 
    
else
    
   error('CHEBFUN:DISKFUNV:rdivide:inputs', 'Unrecognized input arguments.')
   
end
    
end
