function A = useConstraints(disc,blocks)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
L = disc.source;
rows = disc.reproject(blocks);
A = cell2mat(rows);
dim = disc.dimension;
dom = disc.domain;
if ~isempty(L.constraint)
    disc2 = ultraS(L.constraint.operator,dim,dom);
    constr = matrix(disc2);
    A = [ cell2mat(constr); A ];
end
if ~isempty(L.continuity)
    disc2 = ultraS(L.continuity.operator,dim,dom);
    constr = matrix(disc2);
    A = [ cell2mat(constr); A ];
end
end
