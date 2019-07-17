function F = coeffs2spherefunv( X, Y, Z )
%COEFFS2SPHEREFUNV   Convert componentwise matrices of 2D Fourier 
%   coefficients to a spherefunv. 
% 
%   F = coeffs2spherefunv( X, Y, Z ) returns a spherefunv object F that has 
%   matrices of 2D Fourier coefficients X, Y, and Z for each component.  
% 
% See also SPHEREFUN/COEFFS2SPHEREFUN 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

U = spherefun.coeffs2spherefun(X); 
V = spherefun.coeffs2spherefun(Y); 
W = spherefun.coeffs2spherefun(Z); 
F = spherefunv(U, V, W); 
end