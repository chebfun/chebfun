function F = uminus(F)
%UMINUS   Unary minus for a CHEBFUN3T object.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

F.coeffs = -F.coeffs;

end