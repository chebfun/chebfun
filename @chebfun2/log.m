function f = log( f )
%LOG   Natural logarithm of a CHEBFUN2.
%   LOG(F) is the natural logarithm of F. This function does not work if the
%   function passes through or becomes numerically close to zero.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return 
end 

% Positive/negative test. 
[bol, wzero] = singleSignTest( f ); 

if ( bol == 0 ) || ( wzero == 1 )
    error('CHEBFUN:CHEBFUN2:log:notSmooth', ...
    'A change of sign/zero has been detected, unable to represent the result.'); 
end

op = @(x,y) log( feval(f, x, y) );  % Resample.
f = chebfun2( op, f.domain );       % Call constructor. 

end
