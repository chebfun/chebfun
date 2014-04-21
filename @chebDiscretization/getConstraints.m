function B = getConstraints(disc, blocks)
%APPLYCONSTRAINTS   Modify discrete operator to accommodate constraints.
%   M = APPLYCONSTRAINTS(DISC, BLOCKS) uses a cell of matrix BLOCKS created by
%   the discretization DISC in order to return a matrix that incorporates both
%   the linear operator and its constraints.
%
%   [M, P, B, A] = APPLYCONSTRAINTS(DISC, BLOCKS) also returns the individual
%   matrices that represent down-projection of the operator (P), boundary,
%   continuity, and other constraints (B), and the original operator (A). The
%   output M is equivalent to M = [B ; P*A] and is square.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Some preparation.
dim = disc.dimension;
dom = disc.domain;
inputDims = disc.inputDimension(1,:);
discType = str2func( class(disc) );  % constructor of the discretization's type
L = disc.source;

% Now apply the constraints and continuity conditions and collect as rows of B.
B = [];
if ( ~isempty(L.constraint) )
    % Instantiate a discretization of this constraint. 
    disc2 = discType(L.constraint.functional, dim, dom);
    disc2.inputDimension = repmat(inputDims, size(disc2.source, 1), 1);
    constr = matrix(disc2, dim, dom);
    B = [ constr; B ];
end
if ( ~isempty(L.continuity) )
    % Instantiate a discretization of this constraint. 
    disc2 = discType(L.continuity.functional, dim, dom);
    disc2.inputDimension = repmat(inputDims, size(disc2.source, 1), 1);
    constr = matrix(disc2, dim, dom);
    B = [ constr; B ];
end

end
