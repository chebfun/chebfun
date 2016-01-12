function projOrder = getProjOrder(L)
%GETPROJORDER   Get projection order of a discretization.
%   Each boundary and continuity constraint in a linop forces a reduction in the
%   total number of rows in the discrete operator, so that the composite is
%   square. The reduction is found by down-projection of the result of applying
%   the operator.
%
%   PROJORDER = GETPROJORDER(DISC) returns a matrix of dimensions by which each
%   column of the DISC.SOURCE should be down-projected in order to end with a
%   square system. This value is obtained from calling the GETPROJORDER() method
%   of DISC.SOURCE. If no such method exists, PROJORDER = 0 is returned.
%   
%   PROJORDER = DISC.GETPROJORDER(SOURCE) is an equivalent calling sequence.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(L, 'opDiscretization' ) )
    L = L.source;
end

% LINOP objects have nontrivial projection orders. In general, a CHEBMATRIX will
% have a zero projOrder.
if ( isa(L, 'linop') )
    projOrder = getProjOrder(L);
else
    projOrder = 0;
end
    
end
