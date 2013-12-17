function [xx, yy] = chebpts2(nx, ny, varargin)
%CHEBPTS2 Chebyshev tensor points
%
% [XX YY] = CHEBPTS2(N) constructs an N by N grid of Chebyshev 
% tensor points on [-1 1]^2. 
%
% [XX YY] = CHEBPTS2(NX,NY) constructs an NX by NY grid of
% Chebyshev tensor points on [-1 1]^2.  
%
% [XX YY] = CHEBPTS2(NX,NY,D) constructs an NX by NY grid of 
% Chebyshev tensor points on the rectangle [a b] x [c d], where 
% D = [a b c d].
%
% See also CHEBPTS.

if ( nargin > 2 )  % third argument should be a domain. 
   D = varargin{1}; 
   D = D(:).';  % make a row vector.   
   if ( ~all( size( D ) == [1 4] ) )
        error('CHEBFUN2:CHEBPTS2:DOMAIN','Unrecognised domain');
   end
else  % default to the canoncial domain.  
    D = [-1, 1, -1, 1];
end

if ( nargin == 1 ) % make it a square chebyshev grid. 
    ny = nx; 
end

x = chebpts( nx, D(1:2) ); y = chebpts( ny, D(3:4) ); 
[xx, yy] = meshgrid(x, y);   % tensor product. 

end 