function F = diff(F, n, dim)
%DIFF   Componentwise derivative of a CHEBFUN3V object.
%   DIFF(F) is the derivative of each component of F along the x direction.
%
%   DIFF(F, N) is the Nth derivative of each component of F in the x direction.
%
%   DIFF(F, N, DIM) is the Nth derivative of F along the dimension DIM.
%     DIM = 1 (default) is the derivative in the x-direction.
%     DIM = 2 is the derivative in the y-direction.
%     DIM = 3 is the derivative in the z-direction.
%
%   DIFF(F,[NX NY NZ]) is the NX-th partial derivative of F in the first 
%   variable, NY-th derivative of F in the second variable and NZ-th 
%   derivative of F is the third variable. For example, DIFF(F,[2 3 1]) is
%   d^6F/d^2xd^3ydz.

%   DIFF(F,[NX NY]) is the partial derivative of NX of F in the first variable,
%   and NY of F in the second derivative. For example, DIFF(F,[1 2]) is
%   d^3F/dxd^2y.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) ) 
    return
end

% Defaults:
if ( ( nargin == 1 ) || isempty(n) )
    n = 1;
end
if ( nargin < 3 ) 
    dim = 1; 
end

% Diff each component. 
for j = 1:F.nComponents
    F.components{j} = diff(F.components{j}, n, dim); 
end

end