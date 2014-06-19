function outputSpace = getOutputSpace(source)
%GETOUTPUTSPACE    Obtain the range of the ultrapspherical spectral operator.
%   OUT = GETOUTPUTSPACE(SOURCE) returns OUT to indiciate which ultrapsherical
%   polynomial basis to represent the range of the operator. The ultrapsherical
%   polynomial basis is choice so that the differential operator remains sparse.
%   Roughly, if the leading term is of diffOrder K, then the range should be
%   represented in the C^{(K)} polynomial basis.

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
