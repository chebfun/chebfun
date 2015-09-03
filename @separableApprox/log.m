function f = log( f )
%LOG   Natural logarithm of a SEPARABLEAPPROX.
%   LOG(F) is the natural logarithm of F. This function does not work if the
%   function passes through or becomes numerically close to zero.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return 
end 

% Positive/negative test. 
[bol, wzero] = singleSignTest( f ); 

if ( bol == 0 ) || ( wzero == 1 )
    error('CHEBFUN:SEPARABLEAPPROX:log:notSmooth', ...
    'A change of sign/zero has been detected, unable to represent the result.'); 
end

f = compose( f, @log ); 

end
