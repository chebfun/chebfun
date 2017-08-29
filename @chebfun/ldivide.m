function h = ldivide(f,g)
%.\   Pointwise CHEBFUN left divide.
%   F.\G returns a CHEBFUN that represents the function G(x)/F(x). 
%
% See also RDIVIDE, MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call RDIVIDE():
h = rdivide(g, f);

end
