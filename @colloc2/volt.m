function V = volt(disc,kernel,onevar)
%VOLT      Volterra integral operator.
%
%   See also OPERATORBLOCK.VOLT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default onevar to false.
if ( nargin==2 )
    onevar=false; 
end    

% Make use of the cumsum operator. Note that while C(n) would be triangular
% for low-order quadrature, for spectral methods it is not.
C = chebmatrix( operatorBlock.cumsum(disc.domain) );

% Matrix form. Each row of the result, when taken as an inner product with
% function values, does the proper quadrature.
x = points(disc,2);
n = disc.dimension;

if ( onevar )
    V = kernel(x);
else
    [X,Y] = ndgrid(x);
    V = kernel(X,Y);
end

V = V .* matrix(C,n); 

end

