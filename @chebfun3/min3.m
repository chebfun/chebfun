function [val, loc] = min3( f )
%MAX3   Global minimum of a CHEBFUN3.
%
%   See also CHEBFUN3/minandmax3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX3:
[val, loc] = minandmax3(f);

% Extract maximum:
val = val(1);
loc = loc(1, :);

end
