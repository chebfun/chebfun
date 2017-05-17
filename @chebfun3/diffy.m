function g = diffy(f, n)
%DIFFY   Differentiate a CHEBFUN3 with respect to its second argument.
%   G = DIFFY(F) returns a CHEBFUN3 representing the derivative of F in its 
%   second argument. This is the same as DIFF(F, 1, 2).
%
%   G = DIFFY(F, N) returns a CHEBFUN3 representing the Nth derivative of F
%   in its second argument. This is the same as DIFF(F, N, 2).
% 
% See also CHEBFUN3/DIFFX, CHEBFUN3/DIFFZ, and CHEBFUN3/DIFF. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 ) 
    n = 1; 
end

% Call diff:
g = diff(f, n, 2);

end