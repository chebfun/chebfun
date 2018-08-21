function g = spherefun2ballfun(f)
%SPHEREFUN2BALLFUN  Construct BALLFUN whose boundary is prescribed.
%   SPHEREFUN2BALLFUN(F) returns the BALLFUN function whose boundary
%   is equal to F.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

X = coeffs2(f).';
[n,p] = size(X);
G = zeros(2,n,p);
G(1,:,:) = X;
g = ballfun(G);
end
