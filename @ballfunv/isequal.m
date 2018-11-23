function b=isequal(f, g)
%ISEQUAL Equality test for BALLFUNV.  
%   BOL = ISEQUAL(F, G) returns 0 or 1. If returns 1 then F and G are the same
%   BALLFUNV, up to relative machine precision. If returns 0 then F and G are
%   not the same up to relative machine precision. 
%
%   See also ISZERO.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = iszero(f-g);
end
