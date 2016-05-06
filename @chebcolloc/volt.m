function V = volt(kernel, disc, oneVar)
%VOLT    Volterra integral operator for CHEBCOLLOC.
%   For the calling sequence to this method, see OPERATORBLOCK.VOLT.
% 
% See also OPERATORBLOCK.VOLT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default onevar to false.
if ( nargin == 2 )
    oneVar = false; 
end    

% Each row of the result, when taken as an inner product with
% function values, does the proper quadrature.
x = functionPoints(disc);
n = disc.dimension;

if ( oneVar )
    V = kernel(x);
else
    [X, Y] = ndgrid(x);
    V = kernel(X, Y);
end

% Make use of the cumsum operator. Note that while C(n) would be triangular
% for low-order quadrature, for spectral methods it is not.
disc.source = operatorBlock.cumsum(disc.domain);
 
V = V .* matrix(disc, n); 

end

