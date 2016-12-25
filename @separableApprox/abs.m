function f = abs(f)
%ABS   Absolute value of a SEPARABLEAPPROX.
%   ABS(F) returns the absolute value of a SEPARABLEAPPROX.  This function does
%   not work if the function passes through or becomes numerically close to
%   zero.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) ) % check for empty SEPARABLEAPPROX.
    return 
end

if ( isreal(f) )
    % Positive/negative test.
    bool = singleSignTest(f);  % Returns TRUE if there is no sign change.
    if ( ~bool )
        error('CHEBFUN:SEPARABLEAPPROX:abs:notSmooth', ...
            'Sign change detected. Unable to represent the result.'); 
    end
end       

% Still call the constructor in case we missed a change of sign. 
f = compose(f, @abs);

end