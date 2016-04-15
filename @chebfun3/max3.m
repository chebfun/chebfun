function [val, loc] = max3(f)
%MAX3   Global maximum of a CHEBFUN3.
%   See also minandmax3.

% Call MINANDMAX3:
[val, loc] = minandmax3(f);   

% Extract maximum:
val = val(2);
loc = loc(2, :);

end
