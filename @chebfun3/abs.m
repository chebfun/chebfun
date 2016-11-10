function f = abs(f)
%ABS   Absolute value of a CHEBFUN3.
%   ABS(F) returns the absolute value of a CHEBFUN3 object F. This function
%   gives an error if the function passes through or becomes numerically 
%   close to zero.
%
% See also CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) ) % check for empty CHEBFUN3.
    return 
end

if ( isreal(f) )
    % Positive/negative test.
    [ss, ~, ispos] = singleSignTest(f); % Returns TRUE if there is no sign change.
    if ( ss )
        if (ispos)
            return
        else
            f = uminus(f);
            return
        end
    elseif( ~ss )
        error('CHEBFUN:CHEBFUN3:abs:notSmooth', ...
            'Sign change detected. Unable to represent the result.');
    end
else
    % Still call the constructor in case we missed a change of sign. 
    f = compose(f, @abs);
end

end