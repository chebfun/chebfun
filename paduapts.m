function [xy, idx] = paduapts(n, dom)
%PADUAPTS   Padua points.
%   XY = PADUAPTS(N) returns a two-column vector containing the x and y
%   co-ordinates of the degree N first-kind Padua points on [-1 1] x [-1,1].
%
%   XY = PADUAPTS(N, [a, b, c, d]) returns the degree N Padua points on the
%   domain [a b] x [c d].
%
%   The ordering and defniiton of the points is consistent with Padua2DM [1].
%
%   [XY, IDX] = PADUAPTS(...) returns also a logical matrix IDX which denotes
%   the entries of the Chebyshev (N+1)x(N+2) tensor product grid that form XY.
%
% See also CHEBPTS.

% References:
%   [1]: Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello
%        "Padua2DM: fast interpolation and cubature at the Padua points in
%         Matlab/Octave."

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default domain:
if ( nargin < 2 )
    dom = [-1, 1, -1, 1];    
end

% Trivial case:
if ( n == 0 )
    xy = dom([1, 3]);
    return
end

% 1D Chebyshev grids:
xn1 = chebpts(n+1, dom(1:2), 2);
xn2 = chebpts(n+2, dom(3:4), 2);

% Order is flipped from Chebfun for consistency with [1]:
xn1 = xn1(end:-1:1);
xn2 = xn2(end:-1:1);

% Full tensor grid.
[x1, x2] = meshgrid(xn1, xn2); 

% Extract every other term.
idx = true(1, (n+1)*(n+2)); 
idx(1:2:end) = false; 
if ( mod(n, 2) )
    idx = reshape(idx, n+2, n+1);
else
    idx = reshape(idx, n+1, n+2).';
end

% Padua points.
xy = [x1(idx), x2(idx)]; 

end
