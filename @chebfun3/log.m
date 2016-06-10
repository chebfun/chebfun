function f = log(f)
%LOG   Natural logarithm of a CHEBFUN3.
%   LOG(F) is the natural logarithm of F. This function gives an error 
%   message if the function passes through or becomes numerically close to 
%   zero.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    return 
end 

% Positive/negative test. 
[bool, wzero] = singleSignTest(f); 

if ( (bool == 0) || (wzero == 1) )
    error('CHEBFUN:CHEBFUN3:log:notSmooth', ...
    'A change of sign/zero has been detected, unable to represent the result.'); 
end

f = compose(f, @log);

end