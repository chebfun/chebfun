function out = chebcoeffs(f, N)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a TRIGTECH.
%   A = CHEBCOEFFS(F, N) returns the column vector of coefficients such that F
%   = A(1) T_0(x) + ... + A(N-1) T_{N-2}(x) + A(N) T_{N-1}(x), where T_k(x)
%   denotes the k-th Chebyshev polynomial.
%
%   If F is array-valued with P columns, then A is an NxP matrix.
%
% See also LEGCOEFFS, TRIGCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case:
if ( N <= 0 )
    out = [];
    return
end

if ( (nargin == 1) || isempty(N) )
    error('CHEBFUN:TRIGTECH:chebcoeffs:input', ...
        'F does not have a finite Chebyshev series. Please input N.');
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
