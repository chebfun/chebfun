function f = div( f )
%DIV   Divergence of a CHEBFUN3V.
%
%   DIV(F) returns the divergence of the CHEBFUN3V i.e.,
%       divergence(F) = F_x + F_y + F_z.
%
%  This is shorthand for the command DIVERGENCE. 
% 
% See also DIVERGENCE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = divergence(f);

end