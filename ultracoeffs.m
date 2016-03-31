function c = ultracoeffs(f, n, lam)
%ULTRACOEFFS   Compute Ultraspherical series coefficients of a CHEBFUN.
%   A = ULTRACOEFFS(F, N, LAM) returns the first N+1 coefficients in the
%   ultraspherical series expansion of the CHEBFUN F, so that such that F
%   approximately equals A(1) C_0(x) + ... + A(N+1) C_N(x) where C_N(x) denotes
%   the N-th ultraspherical polynomial with parameter LAM. A is a column vector.
%
%   If F is smooth (i.e., numel(f.funs) == 1), then A = ULTRACOEFFS(F, LAM) will
%   assume that N = length(F).
%
%   ULTRACOEFFS does not support quasimatrices.
%
% See also CHEBCOEFFS, JACCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  This is essentially just a wrapper for JACCOEFFS with the correct scaling
%  applied as a post-processing step. Since the recurrence relations and
%  asymptotics for ultraspherical polynomials are easier there would be some
%  benefit from writing a pure ULTRACOEFFS, but for now this should suffice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Deal with Legendre and Chebyshev as special cases.

if ( nargin == 2 )
    lam = n;
    ab = lam - .5;
    c = jaccoeffs(f, ab, ab);
else
    ab = lam - .5;
    c = jaccoeffs(f, n, ab, ab);   
end

n = length(c) - 1;
nn = (0:n).';
% scl = feval(jacpoly(nn, lam-.5, lam-.5), 1)./feval(ultrapoly(nn, lam), 1);
scl = ( gamma(2*lam) ./ gamma(lam+.5) ) * ...
            exp( gammaln(lam+.5+nn) - gammaln(2*lam+nn) );
c = scl.*c;

end
