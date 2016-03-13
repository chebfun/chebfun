function H = ldivide( F, G )
%.\   Pointwise DISKFUNV left divide.
%
% See also RDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(G) ) )
    H = diskfunv;
    return
end

if ( isa( F, 'double' ) )                % DOUBLE .\ DISKFUNV
    nG = G.nComponents;
    Gc = G.components;
    H = G;
    for j = 1:nG
        H.components{j} = ldivide(F, Gc{j});
    end
elseif ( isa(G, 'double') )              % DISKFUNV .\ DOUBLE 
    nF = F.nComponents;
    Fc = F.nComponents;
    H = F;
    for j = 1:nF
        H.components{j} = ldivide(Fc{j}, G);
    end 
else                                     % DISKFUNV.\DISKFUNV 
    nF = F.nComponents;
    nG = G.nComponents;
    if ( nF ~= nG )
        error('DISKFUN:DISKFUNV:ldivide:dimensionMismatch', ...
            'DISKFUNV do not have the same number of components.')
    end
    H = F;
    for j = 1 : nF
        H.components{j} = ldivide( F.components{j}, G.components{j} );
    end
end

end
