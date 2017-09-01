function G = diffy(F, k)
%DIFFY   Differentiate a CHEBFUN3V with respect to its second argument.
%   DIFFY(F) returns a CHEBFUN3V representing the derivative of F in its 
%   second argument. This is the same as DIFF(F, 1, 2).
%
%   DIFFY(F, K) returns a CHEBFUN3V representing the Kth derivative of F in
%   its second argument. This is the same as DIFF(F, K, 2).
% 
% See also CHEBFUN3V/DIFFX, CHEBFUN3V/DIFFZ and CHEBFUN3V/DIFF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % Default to first derivative. 
    k = 1; 
end

G = diff(F, k, 2);

end