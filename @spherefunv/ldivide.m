function H = ldivide( F, G )
%.\   Pointwise SPHEREFUNV left divide.
%
% See also SPHEREFUNV/RDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) || isempty(G) )
    H = spherefunv;
    return
end

if ( isa( F, 'double' ) )                % DOUBLE .\ SPHEREFUNV
    Gc = G.components;
    H = G;
    for j = 1:3
        H.components{j} = ldivide(F, Gc{j});
    end
elseif ( isa(G, 'double') )              % SPHEREFUNV .\ DOUBLE 
    Fc = F.components;
    H = F;
    for j = 1:3
        H.components{j} = ldivide(Fc{j}, G);
    end 
else                                     % SPHEREFUNV .\ SPHEREFUNV 
    H = F;
    for j = 1 : 3
        H.components{j} = ldivide(F.components{j}, G.components{j});
    end
end

end
