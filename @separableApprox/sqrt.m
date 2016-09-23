function f = sqrt( f )
%SQRT   Square root.
%   SQRT(F) returns the square root of a positive SEPARABLEAPPROX F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    return
end

% Positive/negative test.
[bol, wzero] = singleSignTest( f );

if ( ( bol == 0 ) || ( wzero == 1 ) )
    error('CHEBFUN:SEPARABLEAPPROX:sqrt:notSmooth', ...
        'A change of sign/zero has been detected, unable to represent the result.');
end

f = compose( f, @sqrt ); 

end
