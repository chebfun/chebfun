function legvals = legcoeffs2legvals( legcoeffs )
%LEGCOEFFS2LEGVALS  Convert Legendre coeffs to Legendre values.
%   LEGVALS = LEGCOEFFS2LEGVALS(LEGCOEFFS) converts the columns vector
%   LEGCOEFFS representing Legendre coefficients in the Legendre series
%       F(X) = C_LEG(1)*P0(X) + ... + C_LEG(N)*P{N-1}(X), 
%   to a vector LEGVALS representing values of the expansion at LEGPTS, i.e., 
%       LEGVALS = F(legpts(N));
%
% See also CHEBFUN.DLT, LEGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for chebfun/dlt.
legvals = chebfun.dlt( legcoeffs ); 

end
