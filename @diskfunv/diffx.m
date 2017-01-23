function G = diffx( F, n )
%DIFFX   Differentiate a DISKFUNV with respect to its first argument.
% 	DIFFX(F) returns a DISKFUNV representing the derivative of F in its
%   first argument. This is the same as DIFF(F,1,2).
%
% 	DIFFX(F,N) returns a DISKFUNV representing the Nth derivative of F in
%   its first argument. This is the same as DIFF(F,1,N).
%
%   This command is for convenience as the syntax for DIFF, inherited from 
%   the DIFF command for matrices, can be confusing.
% 
% See also DIFFY, DIFF. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to first derivative: 
if ( nargin == 1 )
    n = 1; 
end

% Call DISKFUNV/DIFF:
G = diff( F, 1, n );

end
