function v_leg = chebcoeffs2legvals(c_cheb)
%CHEBCOEFFS2LEGVALS  Convert Chebyshev coefficients to Legendre values. 
%   V_LEG = CHEBCOEFFS2LEGVALS(C_CHEB) converts the vector C_CHEB of Chebyshev
%   coefficients to a vector V_LEG of Legendre values such that
%       C_CHEB(1)*T0(X_k) + ... + C_CHEB(N)*T{N-1} = V_LEG_k, k = 0, ... N-1, 
%   where X_k are the Legendre nodes returned by LEGPTS().
% 
% See also LEGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for chebfun.ndct.
v_leg = chebfun.ndct(c_cheb);

end
