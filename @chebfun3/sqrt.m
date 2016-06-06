function f = sqrt(f)
%SQRT   Square root of a CHEBFUN3 object.
%   SQRT(F) returns the square root of a positive CHEBFUN3 object F.
%
% See also CHEBFUN3/POWER and CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    return
end

% Positive/negative test.
[bool, wzero] = singleSignTest(f);

if ( ( bool == 0 ) || ( wzero == 1 ) )
    error('CHEBFUN:CHEBFUN3:sqrt:notSmooth', ...
        'A change of sign/zero has been detected, unable to represent the result.');
end

f = compose(f, @sqrt); 

end