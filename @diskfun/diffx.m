function g = diffx(f, n)
%DIFFX   Differentiate a DISKFUN with respect to x.
%
%   G = DIFFX(F) returns a DISKFUN representing the derivative of F in its
%   first argument. This is the same as DIFF(F,1,1).
%
%   G = DIFFX(F,N) returns a DISKFUN representing the Nth derivative of F in
%   its first argument. This is the same as DIFF(F,1, N).
%
%   This command is for convenience as the syntax for DIFF, inherited from the
%   DIFF command for matrices, can be confusing.
% 
% See also DISKFUN/DIFFY, DISKFUN/DIFF. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative:
if ( nargin == 1 ) 
    n = 1;
end

% Call diff:
g = diff(f, 1, n);

end
