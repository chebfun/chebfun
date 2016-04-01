function F = vertcat(F , G)
%VERTCAT   Vertical concatenation of SPHEREFUNV objects.
%   [F ; f] where F is a SPHEREFUNV with two components, and f is a SPHEREFUN or
%   scalar then returns a SPHEREFUNV with three components. The first and second
%   component remain unchanged and the third component is f.
% 
%   [f ; F] where F is a SPHEREFUNV with two components, and f is a SPHEREFUN or
%   scalar then returns a SPHEREFUNV with three components. The first is f, and
%   the second and third are the first and second components of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) || isempty(G) )
    F = spherefunv;
    return
end

if ( isa(G, 'double') ) 
    Fc = F.components;
    dom = Fc{1}.domain;
    G = spherefun(G, Fc{1}.domain); 
elseif ( isa(F, 'double') ) 
    Gc = G.components;
    dom = Gc{1}.domain;
    F = spherefun(F, Gc{1}.domain); 
elseif ( isa(G, 'spherefun') )
    if ( ~domainCheck(F.components{1}, G) ) 
        error('SPHEREFUN:SPHEREFUNV:vertcat:domain', 'Inconsistent domains.')
    end
    dom = G.domain; 
else
    error('SPHEREFUN:SPHEREFUNV:vertcat:notSupported', ...
        'Vertical concatenation of these objects is not supported.')
end

if ( isa(F, 'spherefunv') )
    op = [ F.components, {G} ]; 
else
    op = [ {F}, G.components ]; 
end

F = spherefunv(op, dom); 

end
