function g = diffz(f, k)
%DIFFZ   Differentiate a CHEBFUN3 with respect to its third argument.
%   G = DIFFZ(F) returns a CHEBFUN3 representing the derivative of F in its 
%   second argument. This is the same as DIFF(F, 1, 3).
%
%   G = DIFFZ(F, K) returns a CHEBFUN3 representing the Kth derivative of F
%   in its second argument. This is the same as DIFF(F, K, 3).
% 
% See also CHEBFUN3/DIFFX, CHEBFUN3/DIFFY, and CHEBFUN3/DIFF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 ) 
    k = 1; 
end

% Call diff:
g = diff(f, k, 3);

end