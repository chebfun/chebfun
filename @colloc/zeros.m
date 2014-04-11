function Z = zeros(disc)
%ZEROS     Zero operator in COLLOC.
%
%   See also FUNCTIONALBLOCK.ZEROS.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

n = disc.dimension;
Z = zeros(sum(n));

end
