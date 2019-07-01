function legvals = chebvals2legvals(chebvals, kind)
%CHEBVALS2LEGVALS   Convert Chebyshev values to Legendre values
%   LEGVALS = CHEBVALS2LEGVALS(CHEBVALS), converts a vector of CHEBVALS
%   representing values of a Chebyshev expansion at CHEBPTS to a vector LEGVALS
%   representing the expansion evaluated at LEGPTS.
%
%   CHEBVALS2LEGVALS(CHEBVALS, 1) is similar except that the input corresponds
%   to values at first-kind Chebyshev points (i.e., CHEBPTS(N, 1)).
% 
% See also LEGVALS2CHEBVALS, CHEBPTS, LEGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Default to second-kind points:
if ( nargin == 1 )
    kind = 2;
end
if ( kind == 1 )
    chebcoeffs = chebtech1.vals2coeffs(chebvals);
elseif (kind == 2 )
    chebcoeffs = chebtech2.vals2coeffs(chebvals);
else
    error('CHEBFUN:chebvals2legvals:kind', ...
        'Invalid Chebyshev kind. Must be 1 or 2.')
end
    
legvals = chebfun.ndct(chebcoeffs);

end
