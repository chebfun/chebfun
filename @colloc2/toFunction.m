function f = toFunction(disc, values)
%TOFUNCTION Convert COLLOC2 discretization to a CHEBFUN. 
%
%   In case input is matrix valued, the output is an array-valued CHEBFUN, where
%   each column of the CHEBFUN corresponds to a column of the input
%   discretization.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Break into one cell per interval. 
if ( disc.numIntervals > 1 )
    values = mat2cell(values,disc.dimension);
end

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain);

end
