function f = sqrt(f)
%SQRT   Square root.
%   SQRT(F) returns the square root of a positive CHEBFUN3 F.

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