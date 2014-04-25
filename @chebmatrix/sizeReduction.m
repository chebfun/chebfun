function [reduce, d, dRow, dVar] = sizeReduction(L)
%SIZEREDUCTION   Deduce row down-projection dimensions.
%   Each boundary and continuity constraint in a linop forces a reduction in the
%   total number of rows in the discrete operator, so that the composite is
%   square. The reduction is found by down-projection of the result of applying
%   the operator.
%
%   SIZEREDUCTION(L) returns a vector of dimensions by which each row of the
%   system should be down-projected in order to end with a square system. If the
%   differential orders of the variables give the correct total result, they are
%   used; otherwise the reductions are spread as evenly as possible.
%
%   [REDUCE, D, DROW, DVAR] = SIZEREDUCTION(L) also returns the
%   differential order of each variable, the differential order of each
%   system row (equation), and the differential order of each system
%   variable (column).

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO: Adjust documentation. By default chebmatrices are not reduced.

d = L.diffOrder;
dRow = max(d, [], 2);
dVar = max(d, [], 1);
reduce = 0*dVar;

end
