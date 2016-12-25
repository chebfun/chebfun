function F = diff(F, dim, n)
%DIFF   Componentwise derivative of a SPHEREFUNV.
%   DIFF(F) is the tangetial derivative of each component of F in
%   x-direction.
%
%   DIFF(F, DIM) is the first tangential derivative of F along the
%   dimension DIM.
%     DIM = 1 (default) is the derivative in the x-direction.
%     DIM = 2 is the derivative in the y-direction.
%     DIM = 3 is the derivative in the z-direction.
%
%   DIFF(F, DIM, N) is the Nth tangential derivative of each component of F
%   along the dimension specified.
%
%   See also SPHEREFUNV/CURL and SPHEREFUNV/DIV.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check:
if ( isempty(F) ) 
    return
end

% Defaults:
if ( nargin == 1 || isempty(dim) )
    dim = 1;
end
if ( nargin < 3 ) 
    n = 1; 
end

% Diff each component. 
for j = 1:F.nComponents
    F.components{j} = diff(F.components{j}, dim, n); 
end

end
