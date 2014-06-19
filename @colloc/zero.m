function Z = zero(disc)
%ZERO   Zero functional in COLLOC.
%
% See also ZEROS, FUNCTIONALBLOCK.ZERO.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;
Z = zeros(1, sum(n));

end
