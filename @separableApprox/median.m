function g = median( f, dim )
%MEDIAN       Median value of a SEPARABLEAPPROX.
%   G = MEDIAN(F) returns a CHEBFUN G representing the median of the SEPARABLEAPPROX
%   along the y direction, i.e., G = @(x) median( F ( x, : ) ).
%
%   G = MEDIAN(F, DIM) returns a CHEBFUN G representing the median of F along
%   the direction given by DIM, i.e., y-direction if DIM = 1 and x-direction if
%   DIM = 2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

% Default to along the y direction:
if ( nargin == 1 ) 
    dim = 1;
end
dom = f.domain;

% We do not know how to achieve this in an efficient way so we are just going to
% the do the tensor product median.

sample = 2049;
vals = sample( f, sample, sample );   % sample on tensor grid. 
mX = median( vals, dim );                   % discrete median.
mX = mX(:);              

if ( dim == 1 )
    interval = dom( 1:2 );
else
    interval = dom( 3:4 );
end

g = chebfun(mX, interval);    % Form CHEBFUN. 
g = simplify( g, [], 'globaltol' );            % Simplify.

end
