function out = rval(f)
%CHEBTECH.RVAL  Evaluate a CHEBTECH at x = 1.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% We simply sum the coefficients:
out = sum(f.coeffs, 1);

end
