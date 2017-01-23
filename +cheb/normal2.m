function p = normal2(mu, Sigma, dom)
%NORMAL2   Bivariate normal distribution.
%   NORMAL2(MU, SIGMA, DOMAIN) returns the joint probability density for two
%   variables whose means are in the vector MU and whose covariance is the 2x2
%   positive definite matrix SIGMA. The result is a chebfun2 defined on the
%   (finite) domain DOMAIN.
%
%   If DOMAIN is not supplied, the domain is taken to be large enough to make
%   the density essentially zero at the boundary.
%
%   Example:
%      p = cheb.normal2([1,0], [1 1;1 5]);
%      surf(p)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( any(eig(Sigma) < 0) || norm(Sigma - Sigma') > 0 )
    error('CHEB:NORMAL2:covariance:nonSymPosDef', ...
        'Covariance matrix must be symmetric positive definite.')
end

if nargin < 3
    sig = svd(Sigma);
    dom = [ mu(1)+4*[-sig(1) sig(1)], mu(2)+4*[-sig(1) sig(1)] ];
end

x = chebfun2(@(x,y) x,dom);
y = chebfun2(@(x,y) y,dom);
invSig = inv(Sigma);
x = x - mu(1);
y = y - mu(2);
z = invSig(1,1)*x.^2 + (invSig(1,2)+invSig(2,1))*x.*y + invSig(2,2)*y.^2;
const = 2*pi*sqrt(det(Sigma));

p = exp(-z/2)/const;

end
