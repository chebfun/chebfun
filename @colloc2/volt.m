function V = volt(disc,kernel,onevar)
% VOLT  Volterra integral operator.
% V = VOLT(K,D) constructs a chebop representing the Volterra integral
% operator with kernel K for functions in domain D=[a,b]:
%
%      (V*v)(x) = int( K(x,y) v(y), y=a..x )
%
% The kernel function K(x,y) should be smooth for best results.
%
% K must be defined as a function of two inputs X and Y. These may be
% scalar and vector, or they may be matrices defined by NDGRID to represent
% a tensor product of points in DxD.
%
% VOLT(K,D,'onevar') will avoid calling K with tensor product matrices X
% and Y. Instead, the kernel function K should interpret a call K(x) as
% a vector x defining the tensor product grid. This format allows a
% separable or sparse representation for increased efficiency in
% some cases.
%
% Example:
%
% To solve u(x) + x*int(exp(x-y)*u(y),y=0..x) = f(x):
% d = domain(0,2);
% x = chebfun('x',d);
% V = volt(@(x,y) exp(x-y),d);
% u = (1+diag(x)*V) \ sin(exp(3*x));
%
% See also fred, chebop.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default onevar to false
if ( nargin==2 )
    onevar=false; 
end    

% Make use of the cumsum operator. Note that while C(n) would be triangular
% for low-order quadrature, for spectral methods it is not.
C = linop.cumsum(disc.domain);

% Matrix form. Each row of the result, when taken as an inner product with
% function values, does the proper quadrature.
x = points(disc,2);
n = disc.dimension;

if onevar
    V = kernel(x);
else
    [X,Y] = ndgrid(x);
    V = kernel(X,Y);
end
V = V .* matrix(C,n); 

end

