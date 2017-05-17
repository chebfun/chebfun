function out = lval(f)
%CHEBTECH.LVAL  Evaluate a CHEBTECH at x = -1.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Negate the odd terms:
c = f.coeffs;
c(2:2:end,:) = -c(2:2:end,:);
% Sum the resulting coefficients:
out = sum(c, 1);

end
