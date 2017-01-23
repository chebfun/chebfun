function legcoeffs = chebcoeffs2legcoeffs(chebcoeffs)
%CHEBCOEFFS2LEGCOEFFS  Convert Chebyshev coefficients to Legendre coefficients. 
%   C_LEG = CHEBCOEFFS2LEGCOEFFS(C_CHEB) converts the vector C_CHEB of Chebyshev
%   coefficients to a vector C_LEG of Legendre coefficients such that
%       C_CHEB(1)*T0 + ... + C_CHEB(N)*T{N-1} = ...
%           C_LEG(1)*P0 + ... + C_LEG(N)*P{N-1},
%   where P{k} is the degree k Legendre polynomial normalized so that
%   max(|P{k}|) = 1.
% 
% See also CHEB2LEG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for cheb2leg.
legcoeffs = cheb2leg(chebcoeffs);

end
