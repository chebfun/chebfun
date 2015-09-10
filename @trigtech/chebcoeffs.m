function out = chebcoeffs(f, N, kind)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a TRIGTECH.
%   A = CHEBCOEFFS(F, N) or A = CHEBCOEFFS(F, N, 1) returns a column vector of
%   the first N coefficients in the expansion of F in a series of Chebyshev
%   polynomials of the first kind, ordered starting with the coefficient of
%   T_0(x).
%
%   A = CHEBCOEFFS(F, N, 2) does the same but for a series expansion in
%   Chebyshev polynomials of second kind, ordered starting with the coefficient
%   of U_0(x).
%
%   If F is array-valued with P columns, then A is an NxP matrix.
%
% See also LEGCOEFFS, TRIGCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( (nargin < 2) || isempty(N) )
    error('CHEBFUN:TRIGTECH:chebcoeffs:input', ...
        'F does not have a finite Chebyshev series. Please input N.');
end

% Use first-kind Chebyshev polynomials by default.
if ( (nargin < 3) || isempty(kind) )
    kind = 1;
end

% Trivial empty case:
if ( N <= 0 )
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
out = chebcoeffs(chebtech1(@(x) f.feval(x)), N, kind);
    
end
