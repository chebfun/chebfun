function [X, Y] = coeffs2vals(U,V)
% COEFFS2VALS  componentwise conversion of matrix of bivariate Chebyshev
% coefficients to values.
%
%  X,Y = COEFFS2VALS(U, V) converts matrices U and V of bivariate
%  Chebyshev coefficients to matrices of samples X and Y corresponding 
% to values sampled on a 2D tensor Chebyshev grid. 
%
% See also COEFFS2, VALS2COEFFS

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

%%

X = chebfun2.coeffs2vals(U); 
Y = chebfun2.coeffs2vals(V); 
end



