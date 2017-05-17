function [U,H] = poldec(f)
%POLDEC   Polar decomposition of a CHEBFUN2.
%   [U, H] = POLDEC(F) computes chebfun2 objects U and H such that F = U*H.
%   The domain of U is the same as that of F and it is a partial isometry,
%   which means all its singular values are 1 (a finite number) or 0.
%   The domain of H is the square [a,b]x[a,b] where [a,b] is the 
%   x-domain of F, and H is Hermitian positive semidefinite (its
%   eigenvalues are all positive).
%
% See also SVD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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
