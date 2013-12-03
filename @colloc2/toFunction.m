function f = toFunction(disc,values)
if ( disc.numIntervals > 1 )
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    values = mat2cell(values,disc.dimension);
end
f = chebfun(values, disc.domain);
end
