function chebvals = legvals2chebvals(legvals)
%LEGVALS2CHEBVALS   Convert Legendre values to Chebyshev values.
%   CHEBVALS = LEGVALS2CHEBVALS(LEGVALS), converts the vector of LEGVALS
%   representing values of a polynomial at LEGPTS to a vector CHEBVALS
%   representing the same polynomial evaluated at CHEBPTS.
% 
% See also CHEBVALS2LEGVALS, CHEBPTS, LEGPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

legcoeffs = chebfun.idlt(legvals);
chebvals = legcoeffs2chebvals(legcoeffs);

end
