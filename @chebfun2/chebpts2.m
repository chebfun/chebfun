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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin > 2 )  
   % Third argument should be a domain. 
   D = D(:).';  % make a row vector.   
   if ( ~all( size( D ) == [1 4] ) )
        error('CHEBFUN:CHEBFUN2:chebpts2:domain', 'Unrecognised domain.');
   end
else  % Default to the canoncial domain.  
    D = [-1, 1, -1, 1];
end

if ( nargin == 1 ) 
    % Make it a square Chebyshev grid if only one input. 
    ny = nx; 
end

tech = chebfunpref().tech(); 

x = tech.chebpts(nx); 
y = tech.chebpts(ny); 

% scale to domain: 
x = (D(2) - D(1))/2*(x+1) + D(1); 
y = (D(4) - D(3))/2*(y+1) + D(3); 
[xx, yy] = meshgrid(x, y);   % Tensor product. 

end 
