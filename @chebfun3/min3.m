function [val, pos] = min3(f)
%MIN3   Global minimum of a CHEBFUN3.
%   [VAL, POS] = MIN3(F) returns the global minimum VAL and the position 
%   POS of the global minimum of a CHEBFUN3 object F.
%
% See also CHEBFUN3/MINANDMAX3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX3:
[val, pos] = minandmax3(f);

% Extract minimum:
val = val(1);
pos = pos(1, :);

end