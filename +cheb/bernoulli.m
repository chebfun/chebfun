function B = bernoulli(N)
%BERNOULLI   Bernoulli polynomials as chebfuns.
%   B = BERNOULLI(N) returns a quasimatrix of the first N+1 Bernoulli
%   polynomials on [0,1].
%
%   Example (Bernolli numbers):
%      B = cheb.bernoulli(4);
%      format rat, B(0,:)

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = chebfun('x', [0,1]);

% Initalize the constant 1 as the first Bernoulli polynomial.
B = 0*x + 1;
% Compute the requested polynomials.
for j = 1:N
    B(:,j+1) = j*cumsum(B(:,j));
    B(:,j+1) = B(:,j+1) - sum(B(:,j+1));
end

end
