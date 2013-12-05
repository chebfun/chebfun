function f = toFunction(disc,values)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( disc.numIntervals > 1 )
    values = mat2cell(values,disc.dimension);
end
f = chebfun(values, disc.domain);
end
