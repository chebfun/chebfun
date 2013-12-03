function outputSpace = getOutputSpace(source)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
if ( isa(source, 'chebmatrix') )
    diffOrders = getDiffOrder(source);
else
    diffOrders = source.diffOrder;
end
outputSpace = max(max(diffOrders, [], 2) - 1, -1);
end
