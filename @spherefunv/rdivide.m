function F = rdivide( F, G )
%./   Pointwise SPHEREFUNV right divide.
%
% See also LDIVIDE.

% Empty check: 
if ( isempty(F) || isempty(G) )
   F = spherefunv;
   return 
end

% Componentwise divide. 
if ( isa(F, 'spherefunv') && isa(G, 'spherefunv') )    
    for jj = 1 : 3 
        F.components{jj} = rdivide( F.components{jj}, G.components{jj} ); 
    end
    
elseif ( isa(F, 'spherefunv') && ( isa(G, 'spherefun') || isa(G, 'double') ) ) 
    
    for jj = 1 : 3 
        F.components{jj} = rdivide( F.components{jj}, G ); 
    end
    
elseif  ( isa(G, 'spherefunv') && ( isa(F, 'spherefun') || isa(F, 'double') ) ) 
    
    H = G; 
    for jj = 1 : 3
        H.components{jj} = rdivide( F, H.components{jj} ); 
    end
    F = H; 
    
else
    
   error('SPHEREFUN:SPHEREFUNV:rdivide:inputs','Unrecognized input arguments.')
   
end
    
end
