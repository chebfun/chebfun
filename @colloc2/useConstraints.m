function [A,P,B] = useConstraints(disc,blocks)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

L = disc.source;
[rows,P] = disc.reproject(blocks);
P = blkdiag(P{:});

B = [];
if ~isempty(L.constraint)
    disc.source = L.constraint.operator;
    constr = matrix(disc);
    B = [ cell2mat(constr); B ];
end
if ~isempty(L.continuity)
    disc.source = L.continuity.operator;
    constr = matrix(disc);
    B = [ cell2mat(constr); B ];
end
A = cell2mat(rows);
A = [ B; A ];

end
