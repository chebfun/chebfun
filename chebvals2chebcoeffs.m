function chebcoeffs = chebvals2chebcoeffs(chebvals, kind)
%CHEBVALS2CHEBCOEFFS  Convert Chebyshev values to coefficients.
% 	CHEBCOEFFS = CHEBVALS2CHEBCOEFFS(CHEBVALS), converts the column vector
%   CHEBVALS of values on a second-kind Chebyshev grid (i.e, F(CHEBPTS(N)))
%   to a vector CHEBCOEFFS of the Chebyshev coefficients of the series
%       F(X) = C_CHEB(1)*T0(X) + ... + C_CHEB(N)*T{N-1}(X).
%
% 	CHEBVALS2CHEBCOEFFS(CHEBVALS, 1) is similar, but assumes the entries in
% 	CHEBVALS come from evaluating on a first-kind Chebyshev grid, i.e.,
%   F(CHEBPTS(N,1))).
% 
% See also CHEBTECH2.VALS2COEFFS, CHEBPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Default to second-kind points
if ( nargin == 1 )
    kind = 2;
end

if ( kind == 1 )
    % This command is a wrapper for chebtech2/vals2coeffs.
    chebcoeffs = chebtech1.vals2coeffs(chebvals);
elseif ( kind == 2 )
    % This command is a wrapper for chebtech1/vals2coeffs.
    chebcoeffs = chebtech2.vals2coeffs(chebvals);
else
    error('CHEBFUN:chebvals2chebcoeffs:kind', ...
        'Invalid Chebyshev kind. Must be 1 or 2.')
end

end
