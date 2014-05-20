function out = chebpoly(f, N)
%CHEBPOLY   Chebyshev polynomial coefficients of a CHEBTECH.
%   A = CHEBPOLY(F) returns the row vector of coefficients such that F = A(1)
%   T_{N-1}(x) + ... + A(N-1) T_1(x) + A(N) T_0(x), where T_k(x) denotes the
%   k-th Chebyshev polynomial and LENGTH(F) = N. This is equivalent to GET(F,
%   'COEFFS').
%
%   A = CHEBPOLY(F, N) truncates or pads the vector A so that M coefficients of
%   the CHEBTECH F are returned.
%
%   If F is array-valued with M columns, then A is an MxN matrix.
%
% See also LEGPOLY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    N = length(f);
end

out = f.coeffs;
if ( (nargin > 1) && ~isempty(N) )
    [s1, s2] = size(out);
    % Pad:
    out = [zeros(N - s1, s2) ; out];
    % or truncate:
    out(1:(end-N),:) = [];
end

end
