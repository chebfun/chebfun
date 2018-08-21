function g = sum(f)
% SUM Integration of a BALLFUN function over r
%   SUM(f) is the integration of the BALLFUN function f over r

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m,n,p] = size(f);
F = f.coeffs;

% Coefficients of integration of the Chebyshev polynomials
K = zeros(1,m);
for i = 0:m-1
    if mod(i,4)==0
        K(i+1) = -1/(i^2-1);
    elseif mod(i,4)==1
        K(i+1) = 1/(i+1);
    elseif mod(i,4)==2
        K(i+1) = -1/(i^2-1);
    else
        K(i+1) = -1/(i-1);
    end
end

% Matrix of coefficients of Fourier-Fourier after integration over r
A = zeros(n,p);
for i = 1:m
   A = A+K(i).*reshape(F(i,:,:),n,p);
end

% Return the spherefun function of lambda,theta; coeffs2spherefun takes a
% matrix of coefficients of theta and lambda instead of lambda, theta
g = spherefun.coeffs2spherefun(A.');
end
