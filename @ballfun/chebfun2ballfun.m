function g = chebfun2ballfun(f)
%CHEBFUN2BALLFUN   Convert a CHEBFUN in the r-variable to a BALLFUN.
%   CHEBFUN2BALLFUN(f) is the ballfun function g(r,lam,th) = f(r)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Chebyshev coefficients of f
F = chebcoeffs(f);
m = length(F);
% Construction of the ballfun function
G = zeros(m,1,1);
G(:,1,1) = F;
g = ballfun(G);
end
