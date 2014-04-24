function b = rhs(disc, f)
%RHS      Discretize the right-hand side of a linear system.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Developers note: 
%   The original function(s) are discretized, then reduced in dimension the same
%   way as the operator. The constraints are prepended to the top of the vector.

% Create an ULTRAS object from the CHEBMATRIX F:
fDisc = ultraS(f, disc.dimension, disc.domain);
fDisc.outputSpace = disc.outputSpace;

% Instantiate (discretize) the ULTRAS discretisation.
b = cell2mat(instantiate(fDisc));

% Prepend constraints.
L = disc.source;
if ( ~isempty(L.constraint) )
    b = [ L.constraint.values ; b ];
end
if ( ~isempty(L.continuity) )
    b = [ L.continuity.values ; b ];
end

end
