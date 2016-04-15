function [val, loc] = max3(f)
%MAX3   Global maximum of a CHEBFUN3.
%
%   See also minandmax3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX3:
[val, loc] = minandmax3(f);   

% Extract maximum:
val = val(2);
loc = loc(2, :);

end
