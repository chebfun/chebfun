function g = diffz(f, n)
%DIFFZ   Differentiate a CHEBFUN3 with respect to its third argument.
%   G = DIFFZ(F) returns a CHEBFUN3 representing the derivative of F in its 
%   second argument. This is the same as DIFF(F, 1, 3).
%
%   G = DIFFZ(F, N) returns a CHEBFUN3 representing the Nth derivative of F
%   in its second argument. This is the same as DIFF(F, N, 3).
% 
% See also CHEBFUN3/DIFFX, CHEBFUN3/DIFFY, and CHEBFUN3/DIFF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 ) 
    n = 1; 
end

% Call diff:
g = diff(f, n, 3);

end