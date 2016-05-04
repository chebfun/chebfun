function B = getConstraints(disc)
%GETONSTRAINTS   Get constraints and continuity of a linear operator.
%   B = GETONSTRAINTS(DISC) returns a matrix discretization of the constraints
%   and continuity conditions from the linear operator in DISC.SOURCE. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer note:
%   The continuity conditions go above the constraints. This must be the same
%   when the RHS is formed.

% Some preparation:
dom = disc.domain;
dim = disc.dimension;
dimAdjust = disc.dimAdjust(1,:);
discType = str2func(class(disc));
L = disc.source;
B = [];

% Instantiate a matrix discretization of the constraints.
if ( ~isempty(L.constraint) )
    disc2 = discType(L.constraint.functional, dim, dom);
    numRows = size(disc2.source, 1);
    disc2.dimAdjust = repmat(dimAdjust, numRows, 1);
    B = matrix(disc2, dim, dom);
end

% Instantiate a matrix discretization of the continuity conditions.
if ( ~isempty(L.continuity) )
    disc2 = discType(L.continuity.functional, dim, dom);
    numRows = size(disc2.source, 1);
    disc2.dimAdjust = repmat(dimAdjust, numRows, 1);
    B = [  matrix(disc2, dim, dom) ; B ];
end

end
