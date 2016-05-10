function p = chebpolyvalm(p, A)
%CHEBPOLYVALM   Evaluate polynomial with matrix argument. 
%   Y = CHEBPOLYVALM(P, X), when P is a vector of length N+1 whose elements are
%   the Chebyshev coefficients of a polynomial, is the value of the polynomial
%   evaluated with matrix argument X. X must be a square matrix.
%       Y = P(1)*T_N(X) + P(2)*T_{N-1}(X) + ... + P(N)*T_1(X) + P(N+1)*I
%
%   Warning: The matrix X must have a spectrum close to [-1, 1], and the matrix
%   X should not be too non-normal.
% 
% See also POLYVALM, CHEBPOLYVAL.  

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check inputs:
[n, m] = size(A); 
if ( n ~= m )
   error('CHEBFUN:chebpolyvalm:square', 'Matrix must be square'); 
end
if ( ~isvector(p) )
    error('CHEBFUN:chebpolyvalm:vector', 'P must be a vector.'); 
end

% The input vector p is in Matlab style ordering, flip it:
p = flipud(p(:));

% Initialise:
p(1) = 2*p(1);
n = length(p); 
BOld = zeros(m); 
B = zeros(m); 
BNew = B; 
I = eye(m);

% Clenshaw's method.
for k = n:-1:1     
    BOld = B; 
    B = BNew; 
    BNew = p(k)*I + 2*A*B - BOld;
end

% Correct at end.
p = .5*(BNew - BOld);

end
