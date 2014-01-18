function h = ldivide( f, g )
%.\   Pointwise chebfun2 left array divide.
%
% F.\G if G is a chebfun2 and F is a double this returns (1/F)*G
%
% F.\G if G is a double and F is a chebfun2 this returns G\F, but this
% does not work if F becomes numerically close to zero.
%
% F.\G we do not allow F and G to both be chebfun2 object.
% 
% F.\G is the same as the command ldivide(F,G)

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

h = rdivide( g, f );

end