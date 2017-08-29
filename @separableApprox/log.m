function f = log(f)
%LOG   Natural logarithm of a SEPARABLEAPPROX.
%   LOG(F) is the natural logarithm of F. This function returns an error 
%   if the function passes through or becomes numerically close to zero.

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
        error('CHEBFUN:SEPARABLEAPPROX:log:notSmooth', ...
            'Sign change detected. Unable to represent the result.'); 
    end
end

f = compose(f, @log);

end