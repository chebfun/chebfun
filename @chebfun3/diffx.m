function g = diffx(f, n)
%DIFFX   Differentiate a CHEBFUN3 with respect to its first argument.
%
%   G = DIFFX(F) returns a CHEBFUN3 representing the derivative of F in its
%   first argument. This is the same as DIFF(F,1,1).
%
%   G = DIFFX(F,N) returns a CHEBFUN3 representing the Nth derivative of F in
%   its first argument. This is the same as DIFF(F,N,1).
% 
% See also DIFFY, DIFFZ, and DIFF. 


% Default to first derivative:
if ( nargin == 1 ) 
    n = 1;
end

% Call diff:
g = diff(f, n, 1);

end