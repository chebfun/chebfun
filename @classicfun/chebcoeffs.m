function out = chebcoeffs(f, N)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a CLASSICFUN.
%   A = CHEBCOEFFS(F) returns the row vector of coefficients such that F = A(1)
%   T_0(x) + ... + A(M) T_{M - 1}(x) + A(M+1) T_M(x), where T_M(x) denotes the
%   M-th Chebyshev polynomial.
%
%   A = CHEBCOEFFS(F, N) truncates or pads the vector A so that N coefficients of
%   the CLASSICFUN F are returned. If N is not given, N = LENGTH(F) is used by default.
%
%   If F is array-valued with M columns, then A is an MxN matrix.
%
% See also LEGPOLY FOURCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    N = length(f);
end

% Call CHEBCOEFFS() of the .ONEFUN:
out = chebcoeffs(f.onefun, N);

end
