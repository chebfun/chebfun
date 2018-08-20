function g = spherefun2ballfun(f,S)
% SPHEREFUN2BALLFUN BALLFUN function whose boundary is a SPHEREFUN
%   SPHEREFUN2BALLFUN(f, S) is the BALLFUN function whose boundary
%   is equal to f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = S(2); p = S(3);
X = coeffs2(f,n,p);
G = zeros(S);
G(1,:,:) = X.';
g = ballfun(G);
end
