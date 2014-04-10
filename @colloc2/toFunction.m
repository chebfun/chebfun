function f = toFunction(disc, values)
%TOFUNCTION Convert COLLOC2 discretization to a CHEBFUN. 
%   TOFUNCTION(DISC,VALUES) converts the values of a COLLOC2-discretized
%   function to a CHEBFUN. If DISC.DOMAIN has breakpoints, the input should
%   be a vector having the smooth pieces stacked.
%
%   If VALUES is matrix valued, the output is an array-valued CHEBFUN,
%   where each column of the CHEBFUN corresponds to a column of the input.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Break into one cell per interval. 
if ( disc.numIntervals > 1 )
    values = mat2cell(values,disc.dimension);
end

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain, 'chebkind', 1);

end
