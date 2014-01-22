function f = toFunction(disc,values)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( disc.numIntervals > 1 )
    values = mat2cell(values,disc.dimension);
end
% TODO: Just observed this in linop/eigs. Why can we e.g. pass values as 65x33
% matrix, but f will be a one column chebfun (Inf x 33)? AB, 22/1/14.
f = chebfun(values, disc.domain);
end
