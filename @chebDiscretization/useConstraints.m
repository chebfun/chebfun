function [A, P, B] = useConstraints(disc,blocks)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
L = disc.source;
[rows, P] = disc.reproject(blocks);
A = cell2mat(rows);

P = blkdiag(P{:});

dim = disc.dimension;
dom = disc.domain;
discType = str2func( class(disc) );

B = [];
if ~isempty(L.constraint)
    disc2 = discType(L.constraint.operator,dim,dom);
    constr = matrix(disc2);
    B = [ constr; B ];
end
if ~isempty(L.continuity)
    disc2 = discType(L.continuity.operator,dim,dom);
    constr = matrix(disc2);
    B = [ constr; B ];
end

A = [ B; A ];

end
