function f = abs( f )
%ABS Absolute value of a chebfun2.
% 
% ABS(F) returns the absolute value of a chebfun2. This function does 
% not work if the function passes through or becomes numerically close to
% zero. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) ) % check for empty chebfun2.
    return 
end 

% positive/negative test. 
bol = chebfun2.singleSignTest( f );   % returns 1 if there is no sign change

if ( ~bol )
   error('CHEBFUN2:ABS','Sign change detected, unable to represent the result.'); 
end

% Still call the constructor in case we missed a change of sign. 
op = @(x,y) abs( f.feval(x, y) );  % Resample. 
f = chebfun2( op, f.domain );          % Call constructor. 

end