function v_cheb = legcoeffs2chebvals(c_leg, varargin)
%LEGCOEFFS2CHEBVALS  Convert Legendre coefficients to Chebyshev values. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

c_cheb = leg2cheb(c_leg, varargin{:}); % Convert to Chebyshev coefficients.
v_cheb = chebvals2chebcoeffs(c_cheb);  % Convert to Chebyshev values.

end