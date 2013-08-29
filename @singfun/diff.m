function f = diff(f, k)
%DIFF   Derivative of a SINGFUN.
%   DIFF(F) is the derivative of the SINGFUN F, while DIFF(F, K) is its
%   Kth derivative.
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

% Trivial case of an empty SINGFUN:
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

%% Differentiate the SINGFUN F, k times 
while ( k > 0 )
    % Decrease k
    k = k - 1;
        
    % Apply the product rule to the SINGFUN F: f = g(x) .* (1+x).^a .* (1-x).^b
   
    % Three terms of the derivative:
    % fist term: g'(x) .* (1+x).^a .* (1-x).^b
    s = f;
    s.smoothPart = diff(s.smoothPart);
    
    % second term: a * g(x) .* (1+x).^(a-1) .* (1-x).^b    
    % if the exponent at the left end point is non-zero
    if ( f.exponents(1) )
        t = f;    
        t.smoothPart = t.smoothPart * t.exponents(1);
        t.exponents(1) = t.exponents(1)-1;
        s = s + t;
    end
    
    % third term: -b * g(x) .* (1+x).^a .* (1-x).^(b-1)
    if ( f.exponents(2) )
        u = f;
        u.smoothPart = u.smoothPart * (-u.exponents(2));
        u.exponents(2) = u.exponents(2)-1;
        s = s + u;
    end  
    % s is the computed derivative, copy it in f and return
    f = s;
end

end