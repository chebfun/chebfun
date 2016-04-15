function [val, loc] = min3( f )
%MAX3   Global minimum of a CHEBFUN3.
%   See also minandmax3.

% Call MINANDMAX3:
[val, loc] = minandmax3(f);

% Extract maximum:
val = val(1);
loc = loc(1, :);

end
