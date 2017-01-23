function f = sqrt( f )
%SQRT   Square root.
%   SQRT(F) returns the square root of a positive SEPARABLEAPPROX F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    return
end

if ( isreal(f) )
    % Positive/negative test.
    bool = singleSignTest(f);  % Returns TRUE if there is no sign change.
    if ( ~bool )
        error('CHEBFUN:SEPARABLEAPPROX:sqrt:notSmooth', ...
            'Sign change detected. Unable to represent the result.');
    end
end

f = compose(f, @sqrt);

end