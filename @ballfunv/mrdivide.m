function W = mrdivide(V,c)
%/   BALLFUNV right divide.
%
% F/a divides each component of a BALLFUNV F by the scalar a. 
% 
% Only allowed to divide by scalars. 
% 
% See also MTIMES.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

W = V*(1/c);
end
