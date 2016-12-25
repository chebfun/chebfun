function b = rhs(disc, f)
%RHS   Discretize the right-hand side of a linear system for TRIGSPEC.
%   B = RHS(DISC, F) returns a discrete version, B, of the function (or
%   chebmatrix) F, as defined by the discretization DISC.
%
% See also MATRIX, REDUCE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Create a TRIGSPEC object from the CHEBMATRIX F.
fDisc = trigspec(f, disc.dimension, disc.domain);
fDisc.outputSpace = disc.outputSpace;

% Instantiate (discretize) the TRIGSPEC discretisation.
b = cell2mat(instantiate(fDisc));

end
