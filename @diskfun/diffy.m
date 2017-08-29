function g = diffy(f, n)
%DIFFY   Differentiate a DISKFUN f(x,y) with respect to y.
%
%   G = DIFFY(F) returns a DISKFUN representing the derivative of F in y
%   This is the same as DIFF(F,2,1).
%
%   G = DIFFY(F,N) returns a DISKFUN representing the Nth derivative of F in
%   y. This is the same as DIFF(F,2,N).
%
%   This command is for convenience as the syntax for DIFF, inherited from the
%   DIFF command for matrices, can be confusing.
% 
% See also DISKFUN/DIFFX, DISKFUN/DIFF. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 ) 
    n = 1; 
end

% Call diff:
g = diff(f, 2, n);

end
