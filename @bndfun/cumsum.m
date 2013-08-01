function f = cumsum(f, k)
%CUMSUM   Indefinite integral of a BNDFUN.
%   CUMSUM(F) is the indefinite integral of the BNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, K) will compute the Kth indefinite integral with the constant of
%   integration chosen so that each intermediate integral evaluates to 0 at x=a.
%   Thus CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%%
% Trivial case of an empty BNDFUN:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    % Compute first indefinite intergral by default
    k = 1;
end

% Rescaling factor, (b-a)/2, to the mth power.
rescaleFactorK = (.5*diff(f.domain))^k;

% Compute the CUMSUM of all of f's ONEFUNs, multiply by the rescaling factor,
% and assign to the ONEFUN field of f.
f.onefun = cumsum(f.onefun, k)*rescaleFactorK;

end
