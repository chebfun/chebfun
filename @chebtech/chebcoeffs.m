function out = chebcoeffs(f, N)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a CHEBTECH.
%   A = CHEBCOEFFS(F) returns the row vector of coefficients such that F = A(1)
%   T_{N-1}(x) + ... + A(N-1) T_1(x) + A(N) T_0(x), where T_k(x) denotes the
%   k-th Chebyshev polynomial and LENGTH(F) = N. This is equivalent to GET(F,
%   'COEFFS').
%
%   If F is array-valued with P columns, then A is an PxN matrix.
%
%   A = CHEBCOEFFS(F, M) truncates or pads the vector A so that M coefficients
%   of the CHEBTECH F are returned.
%
%   If F is array-valued with P columns, then A is an PxM matrix.
%
% See also LEGCOEFFS, FOURCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    N = length(f);
end

out = f.coeffs;
if ( (nargin > 1) && ~isempty(N) )
    [s1, s2] = size(out);
    % Pad:
    out = [ out ; zeros(N - s1, s2)];
    % or truncate:
    out = out(1:N,:);
end

end