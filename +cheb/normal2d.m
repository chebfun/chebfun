function p = normal2d(mu,Sigma,domn)
%NORMAL2D   Bivariate normal distribution.
%   NORMAL2D(MU,SIGMA,DOMAIN) returns the joint probability density for two
%   variables whose means are in the vector MU and whose covariance is the
%   2x2 positive definite matrix SIGMA. The result is a chebfun2 defined on
%   the (finite) domain DOMAIN.
%
%   If DOMAIN is not supplied, the domain is taken to be large enough to
%   make the density essentially zero at the boundary.
%
%   Example:
%      p = cheb.normal2d([1,0],[1 1;1 5]);
%      surf(p)

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if any(eig(Sigma) < 0)
    error('Covariance matrix must be positive definite.')
end

if nargin < 3
    sig = svd(Sigma);
    domn = [ mu(1)+4*[-sig(1) sig(1)], mu(2)+4*[-sig(1) sig(1)] ];
end

x = chebfun2(@(x,y) x,domn);
y = chebfun2(@(x,y) y,domn);
invSig = inv(Sigma);
x = x - mu(1);
y = y - mu(2);
z = invSig(1,1)*x.^2 + (invSig(1,2)+invSig(2,1))*x.*y + invSig(2,2)*y.^2;
const = 2*pi*sqrt(det(Sigma));

p = exp(-z/2)/const;

end