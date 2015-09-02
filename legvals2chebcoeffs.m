function c_cheb = legvals2chebcoeffs(v_leg)
%LEGVALS2CHEBCOEFFS  Convert Legendre values to Chebyshev coefficients. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

c_leg = chebfun.idlt(v_leg); % Convert to Legendre coefficients.
c_cheb = leg2cheb(c_leg);    % Convert to Chebyshev coefficients.

end