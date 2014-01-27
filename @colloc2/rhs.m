function b = rhs(disc,f)
%RHS      Discretize the right-hand side of a linear system for COLLOC2.

% The original function(s) are discretized, then reduced in dimension the same
% way as the operator. The constraints are prepended to the top of the vector. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

row = instantiate(disc,f.blocks);
row = reduce(disc,row);

b = cell2mat(row);  % transform into matrix

% Prepend the values of the constraints.
L = disc.source;
if ~isempty(L.constraint)
    b = [ L.constraint.values; b ];
end
if ~isempty(L.continuity)
    b = [ L.continuity.values; b ];
end

end
