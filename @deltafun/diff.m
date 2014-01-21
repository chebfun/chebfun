function f = diff(f, k)
%DIFF   Derivative of a DELTAFUN.
%   DIFF(F) is the derivative of the DELTAFUN F, while DIFF(F, K) is its Kth
%   derivative.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%% Check the inputs:

% Trivial case of an empty DELTAFUN:
if ( isempty(f) )
    return
end

if ( (nargin < 2) || isempty(k) )
    % Order of derivative not passed. Assume 1st derivative by default:
    k = 1;
elseif ( k == 0 )
    % Nothing to do here!
    return
end

%% Differentiate the DELTAFUN F, k times
% Differentiate the classical smooth part:
if ( ~isempty(f.funPart) )
    f.funPart = diff(f.funPart, k);
end

% Differentiate the distributional part. This just amounts to shifting
% the magnitude matrix down by k rows by adding k zero rows at the top. 
deltaMag = f.impulses;
m = size(deltaMag, 2);
f.impulses = [ zeros(k, m); deltaMag;];

%%
% Simplify, just in case:
f = simplify(f);                 

end