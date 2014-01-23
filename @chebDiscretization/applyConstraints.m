function [M, P, B, A] = applyConstraints(disc,blocks)


%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Convert the blocks to the original, unconstrained matrix.
A = cell2mat(blocks);

% Project rows down, and record the projection matrix as well.
[rows, P] = disc.reproject(blocks);
M = cell2mat(rows);
P = blkdiag(P{:});

% Now apply the constraints and continuity conditions and collect as rows of B.
dim = disc.dimension;
dom = disc.domain;
discType = str2func( class(disc) );
L = disc.source;

B = [];
if ~isempty(L.constraint)
    disc2 = discType(L.constraint.functional,dim,dom);
    constr = matrix(disc2);
    B = [ constr; B ];
end
if ~isempty(L.continuity)
    disc2 = discType(L.continuity.functional,dim,dom);
    constr = matrix(disc2);
    B = [ constr; B ];
end

% This should restore squareness to the final matrix.
M = [ B; M ];

end
