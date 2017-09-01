function [xx, yy, zz] = chebpts3(nx, ny, nz, dom, kind)
%CHEBPTS3   3D tensor product Chebyshev grid.
%   [XX YY ZZ] = CHEBPTS3(N) constructs an N by N by N grid of Chebyshev 
%   tensor points on [-1 1]^3.
%
%   [XX YY ZZ] = CHEBPTS3(NX,NY,NZ) constructs an NX by NY by NZ grid of 
%   Chebyshev tensor points on [-1 1]^3.
%
%   [XX YY ZZ] = CHEBPTS3(NX,NY,NZ,DOM) constructs an NX by NY by NZ grid 
%   of Chebyshev tensor points on the cube [a b] x [c d] x [e g], where 
%   DOM = [a b c d e g].
%
%   [XX YY ZZ] = CHEBPTS3(NX,NY,NZ,D,KIND) constructor Chebyshev tensor 
%   grid of the kind KIND. KIND = 2 is default.
% 
% See also CHEBPTS and CHEBPTS2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin > 3 )  
   % Third argument should be a domain. 
   dom = dom(:).';  % make a row vector.   
   if ( ~all(size(dom) == [1 6]) )
        error('CHEBFUN:chebpts3:domain', 'Unrecognised domain.');
   end
else  % Default to the canoncial domain.
    dom = [-1, 1, -1, 1, -1, 1];
end

if ( nargin == 1 ) 
    % Make it a square Chebyshev grid if only one input. 
    ny = nx; 
    nz = nx;
end
if ( nargin < 5 ) 
    % default to Chebyshev points of the 2nd kind. 
    kind = 2; 
end

% Get points: 
x = chebpts(nx, dom(1:2), kind); 
y = chebpts(ny, dom(3:4), kind); 
z = chebpts(nz, dom(5:6), kind); 

% Tensor product. 
[xx, yy, zz] = ndgrid(x, y, z); 

end