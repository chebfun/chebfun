function h = ldivide( f, g )
%.\   Pointwise SEPARABLEAPPROX left array divide.
%   F.\G if G is a SEPARABLEAPPROX and F is a double this returns (1/F)*G
%
%   F.\G if G is a double and F is a SEPARABLEAPPROX this returns G\F, but this does
%   not work if F becomes numerically close to zero.
%
%   F.\G we do not allow F and G to both be SEPARABLEAPPROX objects.
% 
%   F.\G is the same as the command ldivide(F,G)

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = rdivide( g, f );

end
