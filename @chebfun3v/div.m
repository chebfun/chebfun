function f = div(f)
%DIV   Divergence of a CHEBFUN3V object.
%   DIV(F) returns the divergence of the CHEBFUN3V i.e., if F = U i + V j + W k,
%   then divergence(F) = U_x + V_y + W_z.
%
%   This is shorthand for the command DIVERGENCE.
% 
% See also CHEBFUN3V/DIVERGENCE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = divergence(f);

end