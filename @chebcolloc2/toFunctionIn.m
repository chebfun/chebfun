function f = toFunctionIn(disc, values)
%TOFUNCTIONIN   Convert CHEBCOLLOC2 discretization to a CHEBFUN. 
%   TOFUNCTIONIN(DISC, VALUES) converts the _rhs_ values of a
%   CHEBCOLLOC2-discretized function (i.e., those at DISC.FUNCTIONONPOINTS) to a
%   CHEBFUN. If DISC.DOMAIN has breakpoints, the input should be a vector
%   having the smooth pieces stacked.
%
%   If VALUES is matrix valued, the output is an array-valued CHEBFUN, where
%   each column of the CHEBFUN corresponds to a column of the input.
%
% See also TOVALUES.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Break into one cell per interval. 
if ( disc.numIntervals > 1 )
    values = mat2cell(values, disc.dimension);
else
    values = {values};
end

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain, 'chebkind', 2);

end
