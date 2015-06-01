function [projOrder, d, dRow, dVar] = getProjOrder(L)
%GETPROJORDER   Get projection order of a LINOP.
%   Each boundary and continuity constraint in a LINOP forces a reduction in the
%   total number of rows in the discrete operator, so that the composite is
%   square. The reduction is found by down-projection of the result of applying
%   the operator.
%
%   GETPROJORDER(DISC) returns a matrix of dimensions by which each column of
%   the system should be down-projected in order to end with a square system.


% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

d = L.diffOrder;
dRow = max(d, [], 2);
dVar = max(d, [], 1);
projOrder = max(dVar, 0);

end
