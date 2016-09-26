function f = diff(f, k)
%DIFF   Derivative of a DELTAFUN.
%   DIFF(F) is the derivative of the DELTAFUN F, while DIFF(F, K) is its Kth
%   derivative.
%
% See also SUM, CUMSUM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case of an empty DELTAFUN:
if ( isempty(f) )
    return
end

if ( (nargin < 2) || isempty(k) )
    % Assume 1st derivative by default:
    k = 1;
elseif ( k == 0 )
    % Nothing to do here!
    return
end

% Differentiate the funPart k times:
f.funPart = diff(f.funPart, k);

% Differentiate the distributional part k times. This just amounts to shifting
% the magnitude matrix down by k rows by adding k zero rows at the top.
deltaMag = f.deltaMag;
m = size(deltaMag, 2);
f.deltaMag = [zeros(k, m) ; deltaMag];

% Simplify:
f = simplifyDeltas(f);                 

end
