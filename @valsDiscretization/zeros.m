function Z = zeros(disc)
%ZEROS    Zero operator in VALSDISCRETIZATION.
%
% See also ZERO, OPERATORBLOCK.ZEROS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;
Z = zeros(sum(n));

end
