function v_cheb = legcoeffs2chebvals(c_leg, varargin)
%LEGCOEFFS2CHEBVALS  Convert Legendre coefficients to Chebyshev values. 
%   V_CHEB = LEGCOEFFS2CHEBVALS(C_LEG) converts the vector C_LEG of Legendre
%   coefficients to a vector V_CHEB of values at second-kind Chebyshev points.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

c_cheb = leg2cheb(c_leg, varargin{:}); % Convert to Chebyshev coefficients.
v_cheb = chebcoeffs2chebvals(c_cheb);  % Convert to Chebyshev values.

end
