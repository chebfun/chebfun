function c_cheb = legvals2chebcoeffs(v_leg)
%LEGVALS2CHEBCOEFFS  Convert Legendre values to Chebyshev coefficients. 
%   C_CHEB = LEGVALS2CHEBCOEFFS(V_LEG) converts the vector V_LEG of values at
%   Legendre points to a vector C_CHEB of Chebyshev coefficients of the unique
%   polynomial that interpolates those values in those points.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

c_leg = chebfun.idlt(v_leg); % Convert to Legendre coefficients.
c_cheb = leg2cheb(c_leg);    % Convert to Chebyshev coefficients.

end
