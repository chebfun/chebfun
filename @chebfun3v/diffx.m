function G = diffx(F, n)
%DIFFX   Differentiate a CHEBFUN3V with respect to its first argument.
% 	DIFFX(F) returns a CHEBFUN3V representing the derivative of F in its
%   first argument. This is the same as DIFF(F,1,1).
%
% 	DIFFX(F, N) returns a CHEBFUN3V representing the Nth derivative of F in
%   its first argument. This is the same as DIFF(F,N,1).
% 
% See also DIFFY, DIFFZ, and DIFF. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    % Default to first derivative. 
    n = 1; 
end

G = diff( F, n, 1 );

end