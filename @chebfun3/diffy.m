function g = diffy(f, k)
%DIFFY   Differentiate a CHEBFUN3 with respect to its second argument.
%   G = DIFFY(F) returns a CHEBFUN3 representing the derivative of F in its 
%   second argument. This is the same as DIFF(F, 1, 2).
%
%   G = DIFFY(F, K) returns a CHEBFUN3 representing the Kth derivative of F
%   in its second argument. This is the same as DIFF(F, K, 2).
% 
% See also CHEBFUN3/DIFFX, CHEBFUN3/DIFFZ, and CHEBFUN3/DIFF. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 ) 
    k = 1; 
end

% Call diff:
g = diff(f, k, 2);

end