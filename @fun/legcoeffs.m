function b = legcoeffs(f, varargin)
%LEGCOEFFS   Compute Legendre series coefficients of a BNDFUN object.
%   B = LEGCOEFFS(F) returns the Legendre series coefficients of BNDFUN F, so
%   that F = B(N+1)*P_N + ... + B(1)*P_0, where P_k is the kth Legendre
%   polynomial (scaled to the domain of F).
%
%   If F is an array-valued BNDFUN, then a matrix of coefficients is returned so
%   that F(:,k) = B(N+1,k)*P_N + ... + B(1,k)*P_0.
%
% See also CHEBCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = legcoeffs(f.onefun, varargin{:});

end
