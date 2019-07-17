function [ X,Y,Z ] = coeffs2vals( U,V,W )
% COEFFS2VALS  componentwise conversion of 2D Fourier coefficient matrices
% to values.
%
%  X,Y, Z = COEFFS2VALS( U, V, W ) converts matrices U, V and W of 2D
%  Fourier coefficients to matrices of samples X,  Y and Z corresponding to
% values sampled on a doubled up tensor grid of equispaced points on the 
% domain [-pi, pi) x [-pi, pi). 
%
% See also COEFFS2, VALS2COEFFS

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

%%

X = spherefun.coeffs2vals(U); 
Y = spherefun.coeffs2vals(V); 
Z = spherefun.coeffs2vals(W); 
end



