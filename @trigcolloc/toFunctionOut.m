function f = toFunctionOut(disc, values, cutoff)
%TOFUNCTIONOUT   Convert TRIGCOLLOC discretization to a CHEBFUN. 
%   TOFUNCTIONOUT(DISC, VALUES, OUT) converts the values of a 
%   TRIGCOLLOC-discretized function to a CHEBFUN. 
%
% See also TOVALUES, TOFUNCTIONIN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for cutoff:
if ( nargin == 3 ) 
    if ( ~mod(cutoff,2) )
        m = cutoff + 1;
    else
        m = cutoff;
    end
else
    m = inf; 
end

% Adjust size of cutoff if necessary:
if ( length(m) ~= numel(values) )
    m = max(m)*ones(numel(values),1);
end

% Cutoff coefficients:
for k = 1:numel(values)
    coeffs = trigtech.vals2coeffs(values);
    [N,~] = size(coeffs);
    if ( ~mod(N,2) )
        coeffs = [.5*coeffs(1,:);coeffs(2:end,:);.5*coeffs(1:end,:)];
        N = N+1;
    end
    C = ceil(N/2);
    m(k) = min(N,m(k));
    c = ceil(m(k)/2);
    inds = 1:m(k);
    inds = inds+(C-c);
    coeffs = coeffs(inds,:);
    values = trigtech.coeffs2vals(coeffs);
end

% Convert VALUES into a CHEBFUN:
f = chebfun(values, disc.domain, 'periodic');

end
