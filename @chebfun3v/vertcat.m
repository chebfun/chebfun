function F = vertcat( F , G )
%VERTCAT Vertical concatenation of CHEBFUN3V objects.
%   [F ; f] where F is a CHEBFUN3V with two components, and f is a CHEBFUN3 or
%   scalar then returns a CHEBFUN3V with three components.  The first and second
%   component remain unchanged and the third component is f.
% 
%   [f ; F] where F is a CHEBFUN3V with two components, and f is a CHEBFUN3 or
%   scalar then returns a CHEBFUN3V with three components. The first is f, and
%   the second and third are the first and second components of F.

if ( isempty( F ) || isempty( G ) )
    F = chebfun3v;
    return
end

if ( isa(G, 'double') ) 
    Fc = F.components;
    dom = Fc{1}.domain;
    G = chebfun3( G, Fc{1}.domain ); 
elseif ( isa(F, 'double') ) 
    Gc = G.components;
    dom = Gc{1}.domain;
    F = chebfun3( F, Gc{1}.domain ); 
elseif ( isa(G, 'chebfun3') )
    if ( ~domainCheck(F.components{1}, G) ) 
        error('CHEBFUN:CHEBFUN3V:vertcat:domain', 'Inconsistent domains.')
    end
    dom = G.domain; 
else
    error('CHEBFUN:CHEBFUN3V:vertcat:notSupported', ...
        'Vertical concatenation of these objects is not supported.')
end

if ( isa(F, 'chebfun3v') )
    op = [ F.components, {G} ]; 
else
    op = [ {F}, G.components ]; 
end

F = chebfun3v( op, dom ); 

end
