function Z = zeros(disc)
%ZEROS    Zero operator in COLLOC.
%
% See also ZERO, OPERATORBLOCK.ZEROS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;
Z = zeros(sum(n));

end
