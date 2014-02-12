function [xx, yy] = chebpts2(nx, ny, D)
%CHEBPTS2 Chebyshev tensor points
%   [XX YY] = CHEBPTS2(N) constructs an N by N grid of Chebyshev tensor points
%   on [-1 1]^2.
%
%   [XX YY] = CHEBPTS2(NX,NY) constructs an NX by NY grid of Chebyshev tensor
%   points on [-1 1]^2.
%
%   [XX YY] = CHEBPTS2(NX,NY,D) constructs an NX by NY grid of Chebyshev tensor
%   points on the rectangle [a b] x [c d], where D = [a b c d].
%
% See also CHEBPTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin > 2 )  
   % Third argument should be a domain. 
   D = D(:).';  % make a row vector.   
   if ( ~all( size( D ) == [1 4] ) )
        error('CHEBFUN2:CHEBPTS2:DOMAIN', 'Unrecognised domain.');
   end
else  % Default to the canoncial domain.  
    D = [-1, 1, -1, 1];
end

if ( nargin == 1 ) 
    % Make it a square Chebyshev grid if only one input. 
    ny = nx; 
end

x = chebpts( nx, D(1:2) );   % x grid.
y = chebpts( ny, D(3:4) );   % y grid
[xx, yy] = meshgrid(x, y);   % Tensor product. 

end 