function f = abs( f )
%ABS Absolute value of a CHEBFUN2.
%   ABS(F) returns the absolute value of a CHEBFUN2. This function does not work
%   if the function passes through or becomes numerically close to zero.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) ) % check for empty CHEBFUN2.
    return 
end 

% Positive/negative test. 
bol = singleSignTest( f );  % Returns TRUE if there is no sign change.

if ( ~bol )
   error('CHEBFUN:CHEBFUN2:abs:notSmooth', ...
       'Sign change detected. Unable to represent the result.'); 
end

% Still call the constructor in case we missed a change of sign. 
op = @(x, y) abs( f.feval(x, y) );   % Resample. 
f = chebfun2( op, f.domain );        % Call constructor. 

end
