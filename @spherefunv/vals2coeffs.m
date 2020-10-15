function [ U, V, W ] = vals2coeffs( X,Y,Z )
% VALS2COEFFS  componentwise conversion of matrices of values to
% matrices of 2D Fourier coefficients.
%
% U, V, W = VALS2COEFFS( X,Y, Z ) converts matrices X, Y and Z of values 
% sampled from doubly periodic functions on equally-spaced tensor grids of 
% the domain [-pi, pi) x [-pi, pi) to matrices U, V and W, 
% containing 2D Fourier coefficients for the corresponding interpolants.
% 
%
% See also COEFFS2VALS

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 


U = spherefun.vals2coeffs(X); 
V = spherefun.vals2coeffs(Y); 
W = spherefun.vals2coeffs(Z); 
end



