function [val, pos] = max3(f)
%MAX3   Global maximum of a CHEBFUN3.
%   [VAL, POS] = MAX3(F) returns the global maximum VAL and the position 
%   POS of the global maximum of a CHEBFUN3 object F.
%
% See also CHEBFUN3/MINANDMAX3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    val = []; 
    pos = [];
    return
end

% Call MINANDMAX3:
[val, pos] = minandmax3(f);   

% Extract maximum:
val = val(2);
pos = pos(2, :);

end
