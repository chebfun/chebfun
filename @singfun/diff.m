function f = diff(f, k, dim)
%DIFF   Derivative of a SINGFUN.
%   DIFF(F) is the derivative of the SINGFUN F, while DIFF(F, K) is its Kth
%   derivative.
%
%   DIFF(F, K) takes the Kth derivative of F.
%
%   Note that the third argument must be 1, indicating that no support for
%   array-valued F.
%
% See also SUM, CUMSUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% Check the inputs:

% Check the dimension, i.e. the third argument:
if ( (nargin == 3) && (dim ~= 1) )
    warning('CHEBFUN:SINGFUN:diff:noSupport', ...
        'SINGFUN does not support array-valued objects.')
end

if ( (nargin < 2) || isempty(k) )
    % Order of derivative not passed in. Assume 1st derivative by default:
    k = 1;
elseif ( k == 0 )
    % Nothing to do here!
    return
end

% Trivial case of an empty SINGFUN:
if ( isempty(f) )
    return
end

%% Differentiate the SINGFUN F, k times 
while ( k > 0 )
    % Decrease k:
    k = k - 1;
        
    % Apply the product rule to the SINGFUN F: f = g(x) .* (1+x).^a .* (1-x).^b
   
    % Three terms of the derivative:
    % First term: g'(x) .* (1+x).^a .* (1-x).^b
    s = f;
    s.smoothPart = diff(s.smoothPart);
    
    % Second term: a * g(x) .* (1+x).^(a-1) .* (1-x).^b    
    if ( f.exponents(1) )
        % If the exponent at the left end point is non-zero.
        t = f;    
        t.smoothPart = t.smoothPart * t.exponents(1);
        t.exponents(1) = t.exponents(1) - 1;
        s = s + t;
    end
    
    % Third term: -b * g(x) .* (1+x).^a .* (1-x).^(b-1)
    if ( f.exponents(2) )
        % If the exponent at the right end point is non-zero.
        u = f;
        u.smoothPart = u.smoothPart * (-u.exponents(2));
        u.exponents(2) = u.exponents(2) - 1;
        s = s + u;
    end 
    
    % s is the computed derivative. Copy it in f:
    f = s;
end

%% 
% If f has negligible exponents return just the SMOOTHPART:
if ( isa(f, 'singfun') && issmooth(f) )
    f = f.smoothPart;
end

end
