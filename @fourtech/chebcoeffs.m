function out = chebcoeffs(f, N)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a FOURTECH.
%   A = CHEBCOEFFS(F) returns the row vector of coefficients such that F = A(1)
%   T_{N-1}(x) + ... + A(N-1) T_1(x) + A(N) T_0(x), where T_k(x) denotes the
%   k-th Chebyshev polynomial and LENGTH(F) = N. This is equivalent to GET(F,
%   'COEFFS').
%
%   A = CHEBCOEFFS(F, N) truncates or pads the vector A so that M coefficients of
%   the FOURTECH F are returned.
%
%   If F is array-valued with M columns, then A is an MxN matrix.
%
% See also LEGCOEFFS, FOURCOEFFS.

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

% TODO: Is there a fast transfrom from Fourier to Chebyshev?

% Since f is a fourtech it is assumed to be smooth and periodic on [-1,1].
% Computing the chebyshev coefficients via innner products requires working with
% non-periodic, but smooth functions on [-1,1]. The right representation for f
% is then a Chebyshev expansion. We therefore convert f to a (happy) Chebyshev
% interpolant and return the coefficients.
f = four2cheb(f);
out = chebcoeffs(f, N);
    
end
