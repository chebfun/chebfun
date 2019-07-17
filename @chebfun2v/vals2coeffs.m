function [U, V] = vals2coeffs(X,Y)
% VALS2COEFFS  componentwise conversion of matrices of values to
% matrices of bivariate Chebyshev coefficients.
%
% U, V = VALS2COEFFS(X, Y) converts matrices X, Y of values representing 
% samples of a function from a 2D Chebyshev tensor grid to matrices U, V 
% of Cbivariate Chebyshev coefficients for the corresponding interpolants.
% 
%
% See also COEFFS2VALS, COEFFS2CHEBFUN2

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

U = chebfun2.vals2coeffs(X); 
V = chebfun2.vals2coeffs(Y); 
end



