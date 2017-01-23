function f = div( f )
%DIV   Divergence of a CHEBFUN2V.
%   DIV(F) returns the divergence of the CHEBFUN2V i.e.,
%       divergence(F) = F_x + F_y.
%
%  This is shorthand for the command DIVERGENCE. 
% 
% See also DIVERGENCE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Note that divergence of a 3-vector is the same, because the functions are
% of two variables.

f = divergence( f ); 

end
