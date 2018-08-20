function g = uminus(f)
% UMINUS BALLFUN unary minus
%   UMINUS(f) is negation of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

X = f.coeffs;
g = ballfun(-X);

end
