function f = log(f)
%LOG   Natural logarithm of a CHEBFUN3.
%   LOG(F) is the natural logarithm of F. This function gives an error 
%   message if the function passes through or becomes numerically close to 
%   zero.
%
% See also CHEBFUN3/COMPOSE.

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
        error('CHEBFUN:CHEBFUN3:log:notSmooth', ...
            'Sign change detected. Unable to represent the result.'); 
    end
end      

f = compose(f, @log);

end