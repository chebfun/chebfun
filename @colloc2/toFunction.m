function f = toFunction(disc, values, inOut)
%TOFUNCTION   Convert COLLOC2 discretization to a CHEBFUN. 
%   TOFUNCTION(DISC, VALUES) or TOFUNCTION(DISC, VALUES, 'OUT') converts
%   the _solution_ values of a COLLOC2-discretized function (i.e., those at
%   DISC.EQUATIONPOINTS) to a CHEBFUN. If DISC.DOMAIN has breakpoints, the
%   input should be a vector having the smooth pieces stacked.
%
%   TOFUNCTION(DISC, VALUES, 'IN') converts the _rhs_ values of a
%   COLLOC2-discretized function (i.e., those at DISC.FUNCTIONONPOINTS) to a
%   CHEBFUN. If DISC.DOMAIN has breakpoints, the input should be a vector
%   having the smooth pieces stacked.
%
%   If VALUES is matrix valued, the output is an array-valued CHEBFUN, where
%   each column of the CHEBFUN corresponds to a column of the input.

% See also TOVALUES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Break into one cell per interval. 
if ( disc.numIntervals > 1 )
    values = mat2cell(values, disc.dimension);
else
    values = {values};
end

if ( nargin < 3 || strcmpi(inOut, 'out') )
    % Convert to values at 2nd-kind points:
    for k = 1:numel(values)
        coeffs = chebtech1.vals2coeffs(values{k});
        values{k} = chebtech2.coeffs2vals(coeffs);
    end
end

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain, 'chebkind', 2);

end
