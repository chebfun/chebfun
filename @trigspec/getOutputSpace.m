function outputSpace = getOutputSpace(source)
%GETOUTPUTSPACE    Output the range of the spectral operator.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(source, 'chebmatrix') )
    diffOrders = getDiffOrder(source);
else
    diffOrders = source.diffOrder;
end

% The output space is determined by the diffOrder information obtained above.
outputSpace = max(max(diffOrders, [], 2) - 1, -1);

end
