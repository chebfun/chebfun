function H = ldivide( F, G )
%.\   Pointwise CHEBFUN3V left divide.
%
% See also RDIVIDE.

if ( ( isempty(F) ) || ( isempty(G) ) )
    H = chebfun3v;
    return
end

if ( isa( F, 'double' ) )                % DOUBLE .\ CHEBFUN3V
    nG = G.nComponents;
    Gc = G.components;
    H = G;
    for j = 1:nG
        H.components{j} = ldivide(F, Gc{j});
    end
    
elseif ( isa(G, 'double') )              % CHEBFUN3V .\ DOUBLE 
    nF = F.nComponents;
    Fc = F.nComponents;
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
    for j = 1 : nF
        H.components{j} = ldivide( F.components{j}, G.components{j} );
    end
end

end