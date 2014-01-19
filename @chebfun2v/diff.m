function F = diff(F, n, dim)
%DIFF Componentwise derivative of a chebfun2v.
%
% DIFF(F) is the derivative of each component of F along the y direction.
%
% DIFF(F,N) is the Nth derivative of each component of F in the y direction.
%
% DIFF(F,N,DIM) is the Nth derivative of F along the dimension DIM.
%     DIM = 1 (default) is the derivative in the y-direction.
%     DIM = 2 is the derivative in the x-direction.
%
% DIFF(F,[NX NY]) is the partial derivative of NX of F in the first 
% variable, and NY of F in the second derivative. For example, DIFF(F,[1
% 2]) is d^3F/dxd^2y.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(F) ) % check for empty chebfun2.
    return
end

if ( nargin == 1 ) % defaults.
    n = 1;
    dim = 1;
end
if ( nargin == 2 ) 
    dim = 1; 
end

if ( isempty(n) )  % empty n defaults to y-derivative.
    n = 1;
end

% Just diff each component. 
for j = 1:F.nComponents
    F.components{j} = diff(F.components{j}, n, dim); 
end

end