function G = diffy( F, n )
%DIFFY   Differentiate a DISKFUNV with respect to its second argument
%   DIFFY(F) returns a DISKFUNV representing the derivative of F in its
%   second argument. This is the same as DIFF(F,2,1).
%
%   DIFFY(F,N) returns a DISKFUNV representing the Nth derivative of F in 
%   its second argument. This is the same as DIFF(F, 2, N).
%
%   This command is for convenience as the syntax for DIFF, inherited from 
%   the DIFF command for matrices, can be confusing.
% 
% See also DIFFX, DIFF. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % Default to first derivative. 
    n = 1; 
end

% Call DISKFUNV/DIFF:
G = diff( F, 2, n);

end
