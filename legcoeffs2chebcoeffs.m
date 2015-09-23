function chebcoeffs = legcoeffs2chebcoeffs(legcoeffs)
%CHEB2LEG  Convert Legendre coefficients to Chebyshev coefficients. 
%   C_CHEB = LEG2CHEB(C_LEG) converts the vector C_LEG of Legendre coefficients
%   to a vector C_CHEB of Chebyshev coefficients such that 
%       C_CHEB(1)*T0 + ... + C_CHEB(N)*T{N-1} = ...
%           C_LEG(N)*P0 + ... + C_LEG(1)*P{N-1}, 
%   where P{k} is the degree k Legendre polynomial normalized so that max(|P{k}|
%   = 1.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for leg2cheb.
chebcoeffs = leg2cheb(legcoeffs);

end
