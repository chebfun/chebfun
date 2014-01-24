function G = diffx( F, n )
%DIFFX differentiate a chebfun2v with respect to its first argument
%
% G = DIFFX(F) returns a chebfun2v representing the derivative of F in its 
% first argument. This is the same as DIFF(F,1,2).
%
% G = DIFFX(F,N) returns a chebfun2v representing the Nth derivative of F in
% its first argument. This is the same as DIFF(F,N,2).
%
% This command is for convenience as the syntax for DIFF, inherited from
% the DIFF command for matrices, can be confusing.
% 
% See also DIFFY, DIFF. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin == 1 )
    n = 1; % default to first derivative. 
end

G = diff( F, n, 2 );

end