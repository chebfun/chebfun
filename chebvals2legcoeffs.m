function c_leg = chebvals2legcoeffs(v_cheb, varargin)
%CHEBVALS2LEGCOEFFS  Convert Chebyshev values to Legendre coefficients. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

c_cheb = chebvals2chebcoeffs(v_cheb);  % Convert to Chebyshev coefficients.
c_leg = cheb2leg(c_cheb, varargin{:}); % Convert to Legendre coefficients.

end