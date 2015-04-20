function g = diffx(f, n)
%DIFFX   Differentiate a CHEBFUN2 with respect to its first argument.
%
%   G = DIFFX(F) returns a CHEBFUN2 representing the derivative of F in its
%   first argument. This is the same as DIFF(F,1,2).
%
%   G = DIFFX(F,N) returns a CHEBFUN2 representing the Nth derivative of F in
%   its first argument. This is the same as DIFF(F,N,2).
%
%   This command is for convenience as the syntax for DIFF, inherited from the
%   DIFF command for matrices, can be confusing.
% 
% See also DIFFY, DIFF. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative:
if ( nargin == 1 ) 
    n = 1;
end

% Call diff:
g = diff(f, n, 2);

end
