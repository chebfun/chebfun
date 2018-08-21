function g = chebfun2ballfun(f,S)
%CHEBFUN2BALLFUN   Convert a CHEBFUN in the r-variable to a BALLFUN.
%   CHEBFUN2BALLFUN(f) is the ballfun function g(r,lam,th) = f(r)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

m = S(1); n = S(2); p = S(3);
% Approximation of f on the m first chebyshev polynomials
f = chebfun(f,m);
% Chebyshev coefficients of f
F = chebcoeffs(f);
% Construction of the ballfun function
G = zeros(m,n,p);
G(:,floor(n/2)+1,floor(p/2)+1) = F;
g = ballfun(G);
end
