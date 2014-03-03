function V = volt(disc, kernel, oneVar)
%VOLT      Volterra integral operator.
%
%   For the calling sequence to this method, see also OPERATORBLOCK.VOLT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default onevar to false.
if ( nargin==2 )
    oneVar = false; 
end    

% Make use of the cumsum operator. Note that while C(n) would be triangular
% for low-order quadrature, for spectral methods it is not.
C = chebmatrix( operatorBlock.cumsum(disc.domain) );

% Matrix form. Each row of the result, when taken as an inner product with
% function values, does the proper quadrature.
x = functionPoints(disc);
n = disc.dimension;

if ( oneVar )
    V = kernel(x);
else
    [X,Y] = ndgrid(x);
    V = kernel(X, Y);
end

V = V .* matrix(C, n); 

end

