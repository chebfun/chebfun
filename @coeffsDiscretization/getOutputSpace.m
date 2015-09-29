function outputSpace = getOutputSpace(source)
%GETOUTPUTSPACE    Obtain the range of the spectral operator.
%   OUT = GETOUTPUTSPACE(SOURCE) returns OUT, the range of the spectral 
%   operator.

%   This is particularly important for the ULTRAS class for which the 
%   ultrapsherical polynomial basis is choosed so that the differential operator 
%   remains sparse. Roughly, if the leading term is of diffOrder K, then the 
%   range should be represented in the C^{(K)} polynomial basis. For TRIGSPEC,
%   the operator is always sparse.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(source, 'chebmatrix') )
    diffOrders = getDiffOrder(source);
else
    diffOrders = source.diffOrder;
end

% The output space is determined by the diffOrder information obtained above.
outputSpace = max(max(diffOrders, [], 2) - 1, -1);

end
