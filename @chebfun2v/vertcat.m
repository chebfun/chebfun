function F = vertcat( F , G )
%VERTCAT Vertical concatenation of chebfun2v objects.
% 
% [F;f] where F is a chebfun2v with two components, and f is a chebfun2
% or scalar then returns a chebfun2v with three components.  The first
% and second component remain unchanged and the third component is f. 
% 
% [f;F] where F is a chebfun2v with two components, and f is a chebfun2 or
% scalar then returns a chebfun2v with three components. The first is f,
% and the second and third are the first and second components of F. 

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information.

if ( isempty( F ) || isempty( G ) )
    F = chebfun2v;
    return
end

if ( isa(G, 'double') ) 
    Fc = F.components;
    G = chebfun2( G, Fc{1}.domain ); 
elseif ( isa(F, 'double') ) 
    Gc = G.components;
    F = chebfun2( F, Gc{1}.domain ); 
elseif ( isa(F, 'chebfun2') || isa(G, 'chebfun2') )
else
        error('CHEBFUN2V:VERTCAT','Vertical concatenation of these objects is not supported.')
end

if ( isa(F, 'chebfun2v') )
    op = [ F.components, {G} ]; 
else
    op = [ {F}, G.components ]; 
end

F = chebfun2v( op ); 

end