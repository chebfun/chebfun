function f = div( f )
%DIV the divergence of a chebfun2v.
%
% F = DIV(F) returns the divergence of the chebfun2v i.e. 
% 
%  div(F) = F_x + F_y
% 
% See also DIVERGENCE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.  

f = divergence(f); 

end