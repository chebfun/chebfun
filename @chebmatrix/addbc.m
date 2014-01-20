function L = addbc(A,varargin)
%ADDBC   Add boundary/general constraint to a chebmatrix.
% 
%   See also LINOP, LINOP.ADDBC.

% Recast the chebmatrix as a linop, then apply the bc.

L = addbc( linop(A), varargin{:} );

end