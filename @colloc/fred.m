function F = fred(disc, kernel, oneVar)
% FRED  Fredholm integral operator using COLLOC2.
%
%   For the calling sequence to this method, see also OPERATORBLOCK.FRED.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default onevar to false
if ( nargin==2 )
    oneVar=false; 
end    

% At given n, multiply function values by CC quadrature weights, then apply
% kernel as inner products.
[x, s] = functionPoints(disc);
n = disc.dimension;

if ( oneVar )  % experimental
    F = kernel(x)*spdiags(s', 0, sum(n), sum(n));
else
    [X, Y] = ndgrid(x);
    F = kernel(X, Y) * spdiags(s', 0,sum(n), sum(n));
end

end