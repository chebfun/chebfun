function out = chebcoeffs(f, N, kind)
%CHEBCOEFFS   Chebyshev polynomial coefficients of a CHEBTECH.
%   A = CHEBCOEFFS(F) returns the column vector of coefficients such that F =
%   A(1) T_0(x) + ... + A(N-1) T_{N-2}(x) + A(N) T_{N-1}(x), where T_k(x)
%   denotes the k-th Chebyshev polynomial and LENGTH(F) = N. This is equivalent
%   to GET(F, 'COEFFS').
%
%   If F is array-valued with P columns, then A is an NxP matrix.
%
%   A = CHEBCOEFFS(F, M) truncates or pads the vector A so that M coefficients
%   of the CHEBTECH F are returned.
%
%   A = CHEBCOEFFS(F, 2) and A = CHEBCOEFFS(F, M, 2) do the same as above but
%   return the coefficients in the expansion of F in the Chebyshev polynomials
%   of the second kind U_k(x).
%
%   If F is array-valued with P columns, then A is an MxP matrix.
%
% See also LEGCOEFFS, FOURCOEFFS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Return number of coefficients equal to length of f by default.
if ( (nargin < 2) || isempty(N) )
    N = length(f);
end

% Use first-kind Chebyshev polynomials by default.
if ( (nargin < 3) || isempty(kind) )
    kind = 1;
end

if ( kind == 1 )
    out = f.coeffs;

    % We only need to pad/truncate if N was specified.
    if ( nargin > 1 )
        [s1, s2] = size(out);
        out = [out ; zeros(N - s1, s2)];  % Pad.
        out = out(1:N,:);                 % Truncate.
    end
elseif ( kind == 2 )
    % We compute 2nd-kind coefficients by computing the 1st-kind coefficients
    % and using CHEBTCOEFFS2CHEBUCOEFFS.  That function uses a recurrence in
    % which the coefficient of U_n requires the coefficients of T_n and T_{n +
    % 2}, so we compute two extra 1st-kind coefficients if 2nd-kind
    % coefficients have been requested.
    out = chebtech.chebTcoeffs2chebUcoeffs(chebcoeffs(f, N + 2, 1));
    out = out(1:N,:);
else
    error('CHEBFUN:CHEBTECH:chebcoeffs:badKind', ...
        '''kind'' input must be 1 or 2.');
end

end
