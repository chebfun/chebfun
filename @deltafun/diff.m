function f = diff(f, k)
%DIFF   Derivative of a DELTAFUN.
%   DIFF(F) is the derivative of the DELTAFUN F, while DIFF(F, K) is its Kth
%   derivative.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%% Check the inputs:

% Trivial case of an empty SINGFUN:
if ( isempty(f) )
    return
end

if ( (nargin < 2) || isempty(k) )
    % Order of derivative not passed in. Assume 1st derivative by default:
    k = 1;
elseif ( k == 0 )
    % Nothing to do here!
    return
end

%% Differentiate the DELTAFUN F, k times 
f.funPart = diff(f.funPart, k);    % Differentiate the classical function.
f.diffOrder = f.diffOrder + k;     % Differentitat the distributional part.

% Make sure the maximum order of the derivative is within limits.
if ( any(f.diffOrder > deltafun.pref.deltafun.maxDiffOrder) )
    error('CHEBFUN:DELTAFUN:DIFF', 'order of highest derivative too high. You can change this through DELTAFUN prefernces.');
end

end
