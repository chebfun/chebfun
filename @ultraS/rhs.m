function b = rhs(disc, f)
% TODO: Documentation missing. AB, 12/2/14.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

fdisc = ultraS(f, disc.dimension, disc.domain);
fdisc.outputSpace = disc.outputSpace;

row = instantiate(fdisc,f.blocks);
row = reduce(disc,row);

b = cell2mat(row);
L = disc.source;
if ~isempty(L.constraint)
    b = [ L.constraint.values; b ];
end
if ~isempty(L.continuity)
    b = [ L.continuity.values; b ];
end

end
