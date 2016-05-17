function [xx, yy] = chebpts2(nx, ny, D, kind)
%CHEBPTS2   Chebyshev tensor product grid.
%   [XX YY] = CHEBPTS2(N) constructs an N by N grid of Chebyshev tensor points
%   on [-1 1]^2.
%
%   [XX YY] = CHEBPTS2(NX,NY) constructs an NX by NY grid of Chebyshev tensor
%   points on [-1 1]^2.
%
%   [XX YY] = CHEBPTS2(NX,NY,D) constructs an NX by NY grid of Chebyshev tensor
%   points on the rectangle [a b] x [c d], where D = [a b c d].
%
%   [XX YY] = CHEBPTS2(NX,NY,D,KIND) constructor Chebyshev tensor grid of
%   the kind KIND. KIND = 2 is default. 
% 
% See also CHEBPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
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
if ( nargin < 4 )
    % default to Chebyshev point of the 2nd kind. 
    kind = 2; 
end

% Get points: 
x = chebpts(nx, D(1:2), kind); 
y = chebpts(ny, D(3:4), kind); 

% Tensor product. 
[xx, yy] = meshgrid(x, y); 

end 
