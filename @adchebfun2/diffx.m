function f = diffx(f, k)
%DIFFX   Differentiate an ADCHEBFUN2 with respect to its first argument.
%
%   G = DIFFX(F) returns an ADCHEBFUN2 representing the derivative of F in its
%   first argument. This is the same as DIFF(F,1,2).
%
%   G = DIFFX(F,N) returns an ADCHEBFUN2 representing the Nth derivative of F in
%   its first argument. This is the same as DIFF(F,N,2).
% 
% See also DIFFY, DIFF. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 2 )
    k = 1; 
end

f = diff( f, k, 2 );

end