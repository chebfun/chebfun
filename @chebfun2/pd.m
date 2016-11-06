function [U,H] = pd(f)
%PD   Polar Decomposition of a CHEBFUN2.
%   U = EIG(F) returns the unitary polar factor of F. 
%   This is the partial isometry with same x-, and y-space as F. 
%   (Note that strictly speaking, U is not unitary but a partial isometry, 
%   with rank(F) copies of singular values at 1. 
%   (as in the so-called "canonical polar decomposition"). 
%
%   [U,H] = EIG(F) returns the unitary polar factor, along with H, a
%   symmetric positive semidefinite bivariate function 
%   (its eigenvalues are all positive). 
%

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


[u, S, v] = svd(f); % SVD of chebfun2
U = u*v'; 
H = v*S*v';

if ( nargout > 1 )
    varargout = { U, H };
else
    varargout = { U };
end

end