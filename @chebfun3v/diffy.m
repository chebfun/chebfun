function G = diffy( F, n )
%DIFFY   Differentiate a CHEBFUN3V with respect to its second argument
%   DIFFY(F) returns a CHEBFUN3V representing the derivative of F in its first
%   argument. This is the same as DIFF(F,1,2).
%
%   DIFFY(F,N) returns a CHEBFUN3V representing the Nth derivative of F in its
%   second argument. This is the same as DIFF(F,N,2).
% 
% See also DIFFX, DIFFZ, and DIFF. 

if ( nargin == 1 ) 
    % Default to first derivative. 
    n = 1; 
end

G = diff( F, n, 2 );

end