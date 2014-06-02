function L = addbc(A, varargin)
%ADDBC   Add boundary/general constraint to a CHEBMATRIX.
% 
% See also LINOP, LINOP.ADDBC.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Recast the chebmatrix as a linop, then apply the bc.
L = addbc(linop(A), varargin{:});

end