function H = ldivide( F, G )
%.\   Pointwise DISKFUNV left divide.
%   F.\G if G is a DISKFUNV and F is a double this returns (1/F)*G
%
%   F.\G if G is a double and F is a DISKFUNV this returns G\F, but this
%   does not work if F becomes numerically close to zero.
%
%   F.\G is not allowed if both F and G are DISKFUNV objects.
% 
%   F.\G is the same as the command LDIVIDE(F, G)
%
% See also RDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( ( isempty(F) ) || ( isempty(G) ) )
    H = diskfunv();
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