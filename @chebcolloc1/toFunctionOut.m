function f = toFunctionOut(disc, values, cutoff)
%TOFUNCTIONOUT   Convert CHEBCOLLOC1 discretization to a CHEBFUN. 
%   TOFUNCTIONIN(DISC, VALUES, CUTOFF) converts the values of a 
%   CHEBCOLLOC1-discretized function to a CHEBFUN. If CUTOFF is
%   specified the resulting CHEBFUN will have length CUTOFF. If DISC.DOMAIN 
%   has breakpoints, the input should have cell arrrays corresponding 
%   to smooth pieces.
%
%   If VALUES is matrix valued, the output is an array-valued CHEBFUN, where
%   each column of the CHEBFUN corresponds to a column of the input.
%
% See also TOVALUES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Break into one cell per interval. 
if ( disc.numIntervals > 1 )
    values = mat2cell(values, disc.dimension);
else
    values = {values};
end

% Check for cutoff
if ( nargin == 3 ) 
    m = cutoff;
else
    m = inf;
end

% adjust size of cutoff if necessary
if ( length(m) ~= numel(values) )
    m = max(m)*ones(numel(values),1);
end

% Cutoff coefficients
for k = 1:numel(values)
    coeffs = chebtech1.vals2coeffs(values{k});
    coeffs = coeffs(1:min(m(k),size(coeffs,1)),:);
    values{k} = chebtech1.coeffs2vals(coeffs);
end

% Convert the VALUES matrix into a CHEBFUN on the appropriate domain
% (potentially an array-valued CHEBFUN).
f = chebfun(values, disc.domain, 'chebkind', 1);

end
