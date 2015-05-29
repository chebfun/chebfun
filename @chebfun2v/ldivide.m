function H = ldivide( F, G )
%.\   Pointwise CHEBFUN2V left divide.
%
% See also RDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(G) ) )
    H = chebfun2v;
    return
end

if ( isa( F, 'double' ) )                % DOUBLE .\ CHEBFUN2V
    nG = G.nComponents;
    Gc = G.components;
    H = G;
    for j = 1:nG
        H.components{j} = ldivide(F, Gc{j});
    end
elseif ( isa(G, 'double') )              % CHEBFUN2V .\ DOUBLE 
    nF = F.nComponents;
    Fc = F.nComponents;
    H = F;
    for j = 1:nF
        H.components{j} = ldivide(Fc{j}, G);
    end 
else                                     % CHEBFUN2V .\ CHEBFUN2V 
    nF = F.nComponents;
    nG = G.nComponents;
    if ( nF ~= nG )
        error('CHEBFUN:CHEBFUN2V:ldivide:dimensionMismatch', ...
            'CHEBFUN2V do not have the same number of components.')
    end
    H = F;
    for j = 1 : nF
        H.components{j} = ldivide( F.components{j}, G.components{j} );
    end
end

end
