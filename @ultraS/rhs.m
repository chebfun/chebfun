function b = rhs(disc, f)
%RHS   Discretize the right-hand side of a linear system for ULTRAS.
%   B = RHS(DISC, F) returns a discrete version, B,  of the function (or
%   chebmatrix) F, as defined by the discretization DISC.
%
% See also MATRIX, REDUCE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developers note: This method works as follows. 
%   The original function(s) are discretized, then reduced in dimension the same
%   way as the operator. The constraints are prepended to the top of the vector.

% Create an ULTRAS object from the CHEBMATRIX F:
fDisc = ultraS(f, disc.dimension, disc.domain);
fDisc.outputSpace = disc.outputSpace;

% Instantiate (discretize) the ULTRAS discretisation.
b = cell2mat(instantiate(fDisc));

% Developer note:
%   The continuity conditions go above the constraints. See getConstraints().

% Prepend the values of the constraints and continuity conditions.
L = disc.source;
if ( ~isempty(L.constraint) )
    b = [ L.constraint.values ; b ];
end
if ( ~isempty(L.continuity) )
    b = [ L.continuity.values ; b ];
end

end
