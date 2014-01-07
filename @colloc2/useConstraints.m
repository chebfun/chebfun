function [A,P,B] = useConstraints(disc,blocks)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

L = disc.source;
[rows,P] = disc.reproject(blocks);
A = cell2mat(rows);

P = blkdiag(P{:});

B = [];
if ~isempty(L.constraint)
    disc.source = L.constraint.operator;
    constr = matrix(disc);
    B = [ constr; B ];
end
if ~isempty(L.continuity)
    disc.source = L.continuity.operator;
    constr = matrix(disc);
    B = [ constr; B ];
end

A = [ B; A ];

end
