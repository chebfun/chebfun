function G = diffx(F, k)
%DIFFX   Differentiate a CHEBFUN3V with respect to its first argument.
%   DIFFX(F) returns a CHEBFUN3V representing the derivative of F in its
%   first argument. This is the same as DIFF(F, 1, 1).
%
% 	DIFFX(F, K) returns a CHEBFUN3V representing the Kth derivative of F in
%   its first argument. This is the same as DIFF(F, K, 1).
% 
% See also CHEBFUN3V/DIFFY, CHEBFUN3V/DIFFZ and CHEBFUN3V/DIFF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    % Default to first derivative. 
    k = 1; 
end

G = diff(F, k, 1);

end