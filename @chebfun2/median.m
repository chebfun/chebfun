function g = median( f, varargin )
%MEDIAN  Median value of a chebfun2
%
% G = MEDIAN(F) returns a chebfun G representing the median of the
% chebfun2 along the y direction, i.e
%
%   G = @(x) median( F ( x, : ) )
%
% G = MEDIAN(F,DIM) returns a chebfun G representing the median of F
% along the direction given by DIM, i.e. y-direction if DIM=1 and
% x-direction if DIM = 2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

% Default to along the y direction:
if ( nargin == 1 ) 
    dim=1;
end

% Take user defined direction:
if ( nargin == 2 ) 
    dim = varargin{1};
end

sample = 2049;
dom = f.domain;

% We do not know how to achieve this in an efficient way so
% we are just going to the do the tensor product median.

vals = chebpolyval2( f, sample, sample );   % sample on tensor grid. 
mX = median( vals, dim );   % discrete median.
mX = mX(:);              

if ( dim == 1 )
    interval = dom( 1:2 );
else
    interval = dom( 3:4 );
end

g = chebfun(mX, interval);    % form chebfun. 
g = simplify( g );            % simplify

end