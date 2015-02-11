function c = truncCoeffs(f, N)
%TRUNCCOEFFS   Trigonometric least square coefficients of a CLASSICFUN.
%   C = TRUNCCOEFFS(F, N) returns the 'middle' 2*N+1 trigonometric 
%   Fourier coefficients of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Call TRUNCCOEFFS() of the .ONEFUN:
c = truncCoeffs(f.onefun, N);
end
