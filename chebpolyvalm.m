function p = chebpolyvalm(c,A)
% CHEBPOLYVALM Evaluate polynomial with matrix argument. 
%
% Y = CHEBPOLYVALM(P,X), when P is a vector of length N+1 whose elements 
% are the Chebyshev coefficients of a polynomial, is the value of the
% polynomial evaluated with matrix argument X.  X must be a square matrix. 
%
%  Y = P(1)*T_N(X) + P(2)*T_{N-1}(X) + ... + P(N)*T_1(X) + P(N+1)*I
%
% Warning: The matrix X must have a spectrum close to [-1,1], and the
% matrix X should not be too non-normal. 
% 
% See also POLYVALM, CHEBPOLYVAL.  

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

[n m]=size(A); 
if (n ~= m)  % square matrix check.
   error('CHEBPOLYVALM:SQUARE','Matrix must be square'); 
end

% flip the coefficients. 
c = c(end:-1:1); 


c(1) = 2*c(1);
n = length(c); Bold = zeros(m); B=zeros(m); Bnew=B; 
% Clenshaw's method.
for k = n:-1:1     
    Bold=B; B=Bnew; 
    Bnew = c(k)*eye(m) + 2*A*B - Bold;
end
% correct at end.
p = .5*(Bnew - Bold);

end