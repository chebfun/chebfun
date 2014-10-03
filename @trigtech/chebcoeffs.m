function out = chebcoeffs(f, N)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a TRIGTECH.
%   A = CHEBCOEFFS(F) returns the row vector of coefficients such that F = A(1)
%   T_{N-1}(x) + ... + A(N-1) T_1(x) + A(N) T_0(x), where T_k(x) denotes the
%   k-th Chebyshev polynomial and LENGTH(F) = N. 
%
%   If F is array-valued with P columns, then A is an PxN matrix.
%
%   A = CHEBCOEFFS(F, M) truncates or pads the vector A so that M coefficients of
%   the TRIGTECH F are returned.
%
%   If F is array-valued with P columns, then A is an PxM matrix.
%
% See also LEGCOEFFS, TRIGCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    N = length(f);
end

% Trivial empty case:
if ( isempty(N) || N <= 0)
    out = [];
    return
end

% [TODO]: Is there a fast transfrom from TRIGTECH to CHEBTECH?
% Since f is a TRIGTECH it is assumed to be smooth and periodic on [-1,1].
% Computing the chebyshev coefficients via innner products requires working with
% non-periodic, but smooth functions on [-1,1]. The right representation for f
% is then a Chebyshev expansion. We therefore convert f to a (happy) Chebyshev
% interpolant and return the coefficients. As an arbitrary choice we will
% convert f to a chebtech1 and then compute the resulting coefficients.
f = chebtech1(@(x) f.feval(x));
out = chebcoeffs(f, N);
    
end