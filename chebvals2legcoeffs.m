function c_leg = chebvals2legcoeffs(v_cheb, varargin)
%CHEBVALS2LEGCOEFFS  Convert Chebyshev values to Legendre coefficients. 
% 	LEGCOEFFS = CHEBVALS2LEGCOEFFS(CHEBVALS), converts the column vector
% 	CHEBVALS of values on a second-kind Chebyshev grid (i.e, F(CHEBPTS(N))) to
% 	a vector LEGCOEFFS of Legendre coefficients, where the degree k Legendre
% 	polynomial P{k} is normalized so that max(|P{k}|) = 1.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

c_cheb = chebvals2chebcoeffs(v_cheb);  % Convert to Chebyshev coefficients.
c_leg = cheb2leg(c_cheb, varargin{:}); % Convert to Legendre coefficients.

end
