function g = diffy(f, n)
%DIFFY   Differentiate a CHEBFUN2 with respect to its second argument.
%
%   G = DIFFY(F) returns a CHEBFUN2 representing the derivative of F in its 
%   second argument. This is the same as DIFF(F,1,1).
%
%   G = DIFFY(F,N) returns a CHEBFUN2 representing the Nth derivative of F in
%   its second argument. This is the same as DIFF(F,N,1).
%
%   This command is for convenience as the syntax for DIFF, inherited from the
%   DIFF command for matrices, can be confusing.
% 
% See also DIFFX, DIFF. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 ) 
    n = 1; 
end

% Call diff:
g = diff(f, n, 1);

end
