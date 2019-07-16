function [X, Y] = vals2coeffs(U,V)
% VALS2COEFFS  componentwise conversion of matrices of values to
% matrices of  Chebyshev-Fourier coefficients.
%
% U, V = VALS2COEFFS(X, Y) converts matrices X, Y of values representing 
% samples of a function from a tensor Chebyshev-Fourier grid 
% and converts them to matrices U, V of Chebyshev-Fourier coefficients 
% for the corresponding interpolant.
% 
%
% See also COEFFS2VALS

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

%%
% Empty check:


X = diskfun.vals2coeffs(U); 
Y = diskfun.vals2coeffs(V); 
end



