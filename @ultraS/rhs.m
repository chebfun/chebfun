function b = rhs(disc, f)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
if isempty(disc.dimension)
    error('Discretization dimension not given.')
end
fdisc = ultraS(f,disc.dimension,disc.domain);
fdisc.outputSpace = disc.outputSpace;
row = matrix(fdisc);
if ~iscell(row)
    row = {row};
end
row = disc.reproject(row);
b = cell2mat(row);
L = disc.source;
if ~isempty(L.constraint)
    b = [ L.constraint.values; b ];
end
if ~isempty(L.continuity)
    b = [ L.continuity.values; b ];
end
end
