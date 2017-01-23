function f = cumsum3(f)
%CUMSUM3   Triple indefinite integral of a CHEBFUN3.
%   F = CUMSUM3(F) returns the triple indefinite integral of a CHEBFUN3. 
%   That is
%                  z  y  x
%                 /  /  /
%   CUMSUM3(F) = |  |  |   F(x,y,z) dx dy dz
%               /  /  /
%              e  c  a
%
%   where [a,b] x [c,d] x [e,g] is the domain of F.
% 
% See also CHEBFUN3/CUMSUM, CHEBFUN3/CUMSUM2, CHEBFUN3/SUM, CHEBFUN3/SUM2 
% and CHEBFUN3/SUM3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) ) 
    f = [];
    return
end

f.cols = cumsum(f.cols);     % CUMSUM along the cols.
f.rows = cumsum(f.rows);     % CUMSUM along the rows.
f.tubes = cumsum(f.tubes);   % CUMSUM along the tubes.

end