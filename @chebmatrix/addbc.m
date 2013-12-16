function L = addbc(A,varargin)

% Recast the chebmatrix as a linop, then apply the bc.

L = addbc( linop(A), varargin{:} );

end