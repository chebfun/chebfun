function f = uminus(f)
%-   BALLFUN unary minus.
%   -F negates the BALLFUN F.
%
%   G = UMINUS(F) is called for the syntax '-F'.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.coeffs = -f.coeffs;

end
