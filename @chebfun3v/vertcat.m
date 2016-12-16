function F = vertcat(F , G)
%VERTCAT   Vertical concatenation of CHEBFUN3V objects.
%
%   [F; f] returns a CHEBFUN3V with three components, if F is a CHEBFUN3V 
%   with two components, and f is a CHEBFUN3 object or a scalar.
% 
%   [f; F] is formed similarly.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) || isempty(G) )
    F = chebfun3v();
    return
end

if ( isa(G, 'double') ) 
    Fc = F.components;
    dom = Fc{1}.domain;
    G = chebfun3(G, Fc{1}.domain);
    
elseif ( isa(F, 'double') )
    Gc = G.components;
    dom = Gc{1}.domain;
    F = chebfun3(F, Gc{1}.domain);
    
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
    op = [F.components, {G}];
else
    op = [{F}, G.components];
end

F = chebfun3v(op, dom);

end