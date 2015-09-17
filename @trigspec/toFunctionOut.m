function f = toFunctionOut(disc, coeffs, cutoff)
%TOFUNCTIONOUT   Convert TRIGSPEC discretization to a CHEBFUN.
%   TOFUNCTIONOUT(DISC, VALUES, OUT) converts the values of a 
%   TRIGSPECC-discretized function to a CHEBFUN. 
%
% See also TOVALUES, TOFUNCTIONIN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for cutoff:
if ( nargin == 3 ) 
    if ( ~mod(cutoff, 2) )
        m = cutoff + 1;
    else
        m = cutoff;
    end
else
    m = inf; 
end

% Adjust size of cutoff if necessary:
if ( length(m) ~= numel(coeffs) )
    m = max(m)*ones(numel(coeffs), 1);
end

% Cutoff coefficients:
for k = 1:numel(coeffs)
    [N, ~] = size(coeffs);
    if ( ~mod(N, 2) )
        coeffs = [.5*coeffs(1,:); coeffs(2:end,:); .5*coeffs(1:end,:)];
        N = N + 1;
    end
    C = ceil(N/2);
    m(k) = min(N, m(k));
    c = ceil(m(k)/2);
    inds = 1:m(k);
    inds = inds + (C - c);
    coeffs = coeffs(inds, :);
end

% Get the domain:
dom = disc.domain; 

% Create a TRIGTECH object from the coefficients: 
tech = trigtech({[], coeffs});

% Create a BNDFUN from the TRIGTECH:
fun{1} = bndfun(tech, struct('domain', dom));

% Create a CHEBFUN from the BNDFUN:
f = chebfun(fun);

end
