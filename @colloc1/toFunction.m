function f = toFunction(disc, values)
%TOFUNCTION   Convert COLLOC1 discretization to a CHEBFUN. 
%   TOFUNCTION(DISC, VALUES) converts the values of a COLLOC1-discretized
%   function to a CHEBFUN. If DISC.DOMAIN has breakpoints, the input should have
%   cell arrrays corresponding to smooth pieces.
%
%   If VALUES is matrix valued, the output is an array-valued CHEBFUN, where
%   each column of the CHEBFUN corresponds to a column of the input.

% See also TOVALUES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Break into one cell per interval. 
if ( disc.numIntervals > 1 )
    values = mat2cell(values, disc.dimension);
end

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain, 'chebkind', 1);

end
