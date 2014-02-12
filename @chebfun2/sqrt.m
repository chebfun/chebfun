function f = sqrt( f )
%SQRT   Square root.
%   SQRT(F) returns the square root of a positive CHEBFUN2 F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    f = chebfun2();
    return
end

% Positive/negative test.
[bol, wzero] = singleSignTest( f );

if ( ( bol == 0 ) || ( wzero == 1 ) )
    error('CHEBFUN2:SQRT', ...
        'A change of sign/zero has been detected, unable to represent the result.');
end

% Call the constructor:
op = @(x,y) sqrt( feval( f, x, y ) ); % Resample.
f = chebfun2( op, f.domain );         % Call constructor.

end