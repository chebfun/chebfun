function H = ldivide(F, G)
%.\   Pointwise left divide for two CHEBFUN3V objects.
%
% See also CHEBFUN3V/RDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(G) ) )
    H = chebfun3v;
    return
end

if ( isa(F, 'double') || isa(F, 'chebfun3') ) % CHEBFUN3 .\ CHEBFUN3V
    nG = G.nComponents; 
    Gc = G.components;
    H = G;
    for j = 1:nG
        H.components{j} = ldivide(F, Gc{j});
    end
    
elseif ( isa(G, 'double') || isa(G, 'chebfun3'))  % CHEBFUN3V .\ CHEBFUN3
    nF = F.nComponents;
    Fc = F.components;
    H = F;
    for j = 1:nF
        H.components{j} = ldivide(Fc{j}, G);
    end 
    
else                                     % CHEBFUN3V .\ CHEBFUN3V 
    nF = F.nComponents;
    nG = G.nComponents;
    if ( nF ~= nG )
        error('CHEBFUN:CHEBFUN3V:ldivide:dimensionMismatch', ...
            'CHEBFUN3V do not have the same number of components.')
    end
    H = F;
    for j = 1:nF
        H.components{j} = ldivide(F.components{j}, G.components{j});
    end
end

end