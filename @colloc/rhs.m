function b = rhs(disc, f)
%RHS      Discretize the right-hand side of a linear system for COLLOC2.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Developers note: This method works as follows. 
%   The original function(s) are discretized, then reduced in dimension the same
%   way as the operator. The constraints are prepended to the top of the vector.

% % Instantiate the RHS on the grid, then reduce
% row = instantiate(disc, f.blocks);
% row = reduce(disc, row);
% b = cell2mat(row);  % transform into matrix

xOut = equationPoints(disc);
for k = 1:numel(f.blocks)
    if ( ~isnumeric(f.blocks{k}) )
        f.blocks{k} = feval(f.blocks{k}, xOut);
    end
end
b = cell2mat(f.blocks);       

% Prepend the values of the constraints. First do the constraints, then add
% continuity conditions.
L = disc.source;
if ( ~isempty(L.constraint) )
    b = [ L.constraint.values; b ];
end
if ( ~isempty(L.continuity) )
    b = [ L.continuity.values; b ];
end

end
