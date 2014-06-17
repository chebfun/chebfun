function out = jaccoeffs(f, n, a, b)
%LEGPOLY   Jacobi polynomial coefficients of a CHEBFUN.
%   A = JACCOEFFS(F, N, ALPHA, BETA) returns the first N+1 coefficients in the
%   Jacobi series expansion of the CHEBFUN F, so that such that F approximately
%   equals A(1) J_N(x) + ... + A(N) J_1(x) + A(N+1) J_0(x) where J_N(x) denotes
%   the N-th Jacobi polynomial with parameters ALPHA and BETA. A is a row
%   vector.
%
%   If F is smooth (i.e., numel(f.funs) == 1), then A = JACCOEFFS(F, ALPHA,
%   BETA) will assume that N = length(F) - 1;
%
%   There is also a JACCOEFFS command in the Chebfun root directory which
%   computes the CHEBFUN corresponding to the Jacobi polynomial J_n(x, ALPHA,
%   BETA). Both versions of JACCOEFFS use the same normalization.
%
%   JACCOEFFS does not support quasimatrices.
%
% See also CHEBCOEFFS, LEGCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    error('CHEBFUN:CHEBFUN:jaccoeffs:quasi', ...
        'JACCOEFFS does not support quasimatrices.');
end

if ( (numel(f.funs) == 1) && (nargin < 4) )
    b = a;
    a = n;
    n = length(f);
end
    
% Jacobi-Vandermonde matrix:
Enorm = jacpoly(0:n-1, a, b, f.domain);

% Compute the projection (with correct scaling):
out = Enorm \ f;

% Convention is that the output is a row vector:
out = out.';

end
                
