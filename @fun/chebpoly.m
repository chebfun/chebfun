function out = chebpoly(f, N)
%CHEBPOLY   Chebyshev polynomial coefficients of a FUN.
%   A = CHEBPOLY(F) returns the row vector of coefficients such that F = A(1)
%   T_M(x) + ... + A(M) T_1(x) + A(M+1) T_0(x), where T_M(x) denotes the M-th
%   Chebyshev polynomial.
%
%   A = CHEBPOLY(F, N) truncates or pads the vector A so that N coefficients of
%   the FUN F are returned. If N is not given, N = LENGTH(F) is used by default.
%
%   If F is array-valued with M columns, then A is an MxN matrix.
%
% See also LEGPOLY.

if ( nargin == 1 )
    N = length(f);
end

% [TODO]: Remove this. (Ths will eventually be implemented by @CHEBTECH/CHEBPOLY
% but we need to wait until changes from upstream have been merged.)
if ( isa(f.onefun, 'chebtech') )
    out = get(f.onefun, 'coeffs');
    if ( (nargin > 1) && ~isempty(N) )
        % Pad / truncate:
        out = [zeros(1, N - length(out)) out];
        out = out(end-(N-1):end);
    end
    return
end

% Call CHEBPOLY() of the .ONEFUN:
out = chebpoly(f.onefun, N);

end
