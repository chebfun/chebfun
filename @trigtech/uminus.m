function f = uminus(f)
%UMINUS   Unary minus of a TRIGTECH.
%   UMINUS(F) is the negative of F.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Negate the values:
f.values = -f.values; 

% Negate the coefficients:
f.coeffs = -f.coeffs;

end
