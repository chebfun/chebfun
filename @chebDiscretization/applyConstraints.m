function [M, P, B, A] = applyConstraints(disc, blocks)
%APPLYCONSTRAINTS Modify discrete operator to accommodate constraints.
%   M = APPLYCONSTRAINTS(DISC, BLOCKS) uses a cell of matrix BLOCKS created by
%   the discretization DISC in order to return a matrix that incorporates both
%   the linear operator and its constraints.
%
%   [M, P, B, A] = APPLYCONSTRAINTS(DISC, BLOCKS) also returns the individual
%   matrices that represent down-projection of the operator (P), boundary,
%   continuity, and other constraints (B), and the original operator (A). The
%   output M is equivalent to M = [B; P*A] and is square.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Convert the blocks to the original, unconstrained matrix.
A = cell2mat(blocks);

% Project rows down, and record the projection matrix as well.
[rows, P] = disc.reduce(blocks);
M = cell2mat(rows);
P = blkdiag(P{:});

% Some preparation.
dim = disc.dimension;
dom = disc.domain;
discType = str2func( class(disc) );  % constructor of the discretization's type
L = disc.source;

% Now apply the constraints and continuity conditions and collect as rows of B.
B = [];
if ( ~isempty(L.constraint) )
    % Instantiate a discretization of this constraint. 
    disc2 = discType(L.constraint.functional, dim, dom);
    constr = matrix(disc2);
    B = [ constr; B ];
end
if ( ~isempty(L.continuity) )
    % Instantiate a discretization of this constraint. 
    disc2 = discType(L.continuity.functional, dim, dom);
    constr = matrix(disc2);
    B = [ constr; B ];
end

% This should restore squareness to the final matrix.
M = [ B; M ];

end
