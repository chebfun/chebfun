function f = sqrt(f)
%SQRT   Square root of a CHEBFUN3 object.
%   SQRT(F) returns the square root of a positive CHEBFUN3 object F.
%
% See also CHEBFUN3/POWER and CHEBFUN3/COMPOSE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    return
end

if ( isreal(f) )
    % Positive/negative test.
    ss = singleSignTest(f);  % Returns TRUE if there is no sign change.
    if ( ~ss )
        error('CHEBFUN:CHEBFUN3:sqrt:notSmooth', ...
            'Sign change detected. Unable to represent the result.'); 
    end
end      

f = compose(f, @sqrt); 

end