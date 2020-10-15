function F = coeffs2diskfunv( X, Y )
%COEFFS2DISKFUNV   Convert componentwise matrices of Chebyshev-Fourier 
%   coefficients to a diskfunv. 
% 
%   F = coeffs2diskfunv( X, Y ) returns a diskfunv object F that has a
%   Chebyshev-Fourier matrices of coefficients X, Y, for each component.  
%   This is useful for computing quantities on the disk with the function F.
% 
% See also DISKFUN/COEFFS2DISKFUN 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

U = diskfun.coeffs2diskfun(X); 
V = diskfun.coeffs2diskfun(Y); 
F = diskfunv(U, V); 
end