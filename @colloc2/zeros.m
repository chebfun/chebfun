function Z = zeros(disc)
%ZEROS     Zero operator in COLLOC2.
%
%   See also FUNCTIONALBLOCK.ZEROs.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

n = disc.dimension;
Z = zeros(sum(n));

end
