function F = coeffs2chebfun2v(X, Y)
%COEFFS2DISKFUNV   Convert componentwise matrices of bivariate Chebyshev
%   coefficients to a chebfun2v. 
% 
%   F = coeffs2chebfun2v(X, Y) returns a chebfun2v object F that has a
%   bivariate Chebyshev matrices of coefficients X, Y, for each component.  
% 
% See also CHEBFUN2V/COEFFS2.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

U = chebfun2(X, 'coeffs'); 
V = chebfun2(Y, 'coeffs'); 
F = chebfun2v(U, V); 
end