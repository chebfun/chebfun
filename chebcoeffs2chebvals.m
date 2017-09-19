function chebvals = chebcoeffs2chebvals(chebcoeffs, kind)
%CHEBCOEFFS2CHEBVALS  Convert Chebyshev coefficients to values.
%   CHEBVALS = CHEBCOEFFS2CHEBVALS(CHEBCOEFFS), converts the column vector
%   CHEBCOEFFS of Chebyshev coefficients in the Chebyshev series
%       F(X) = C_CHEB(1)*T0(X) + ... + C_CHEB(N)*T{N-1}(X), 
%   to a vector CHEBVALS of values of the expansion at second-kind Chebyshev
%   points, i.e., CHEBVALS = F(CHEBPTS(N));
%
% 	CHEBCOEFFS2CHEBVALS(CHEBVALS, 1) is similar, but returns the entries in
% 	CHEBVALS from evaluating on a first-kind Chebyshev grid, i.e.,
%   CHEBVALS = F(CHEBPTS(N,1))). 
% 
% See also CHEBTECH2.COEFFS2VALS, CHEBPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Default to first kind:
if ( nargin == 1 )
    kind = 2;
end

if ( kind == 1 )
    % This command is a wrapper for chebtech2/coeffs2vals:
    chebvals = chebtech1.coeffs2vals(chebcoeffs);
elseif ( kind == 2 )
    chebvals = chebtech2.coeffs2vals(chebcoeffs);
else
    error('CHEBFUN:chebcoeffs2chebvals:kind', ...
        'Invalid Chebyshev kind. Must be 1 or 2.')
end
    

end
