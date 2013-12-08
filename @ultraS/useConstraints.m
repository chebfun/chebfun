function [A, P, B] = useConstraints(disc,blocks)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
L = disc.source;
[rows, P] = disc.reproject(blocks);
P = blkdiag(P{:});

dim = disc.dimension;
dom = disc.domain;
B = [];
if ~isempty(L.constraint)
    disc2 = ultraS(L.constraint.operator,dim,dom);
    constr = matrix(disc2);
    B = [ cell2mat(constr); B ];
end
if ~isempty(L.continuity)
    disc2 = ultraS(L.continuity.operator,dim,dom);
    constr = matrix(disc2);
    B = [ cell2mat(constr); B ];
end

A = cell2mat(rows);
A = [ B; A ];

end
