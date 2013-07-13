function f = diff(f, k)
%DIFF   Derivative of a SINGFUN.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the Kth derivative.
%   
%   [TODO]: Will this make sense in SINGFUN:
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this
%   is a finite difference along the columns of an array-valued CHEBTECH.
%   If F has L columns, an empty CHEBTECH will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.


%% Check the inputs:

% Trivial case of an empty CHEBTECH:
if ( isempty(f) )
    return
end

if ( nargin < 2 || isempty(k) )
    % Order of derivative not passed in. Assume 1st derivative by default:
    k = 1;
elseif ( k == 0 )
    % Nothing to do here!
    return
end

%% Differentiate k times the SINGFUN F
while ( k > 0 )
    % Decrease k
    k = k - 1;
        
    % Apply the product rule to
    % the SINGFUN F:
    % f = g .* (1-x).^a .* (1+x).^b
   
    % Three terms of the derivative:
    % fist term: g' .* (1-x).^a .* (1+x).^b
    s = f;
    s.smoothPart = diff(s.smoothPart);
    
    % second term: -a * g .* (1-x).^(a-1) .* (1+x).^b
    t = f;
    t.smoothPart = t.smoothPart * (-t.exponents(1));
    t.exponents(1) = t.exponents(1)-1;
    
    % third term: b * g .* (1-x).^a .* (1+x).^(b-1)
    u = f;
    u.smoothPart = u.smoothPart * u.exponents(2);
    u.exponents(2) = u.exponents(2)-1;
    
    % The derivative is the sum of
    % the above three functions.
    f = s + t + u;
end
end