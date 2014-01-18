function f = exp(f) 
% EXP  Exponential of a chebfun2
%
% EXP(F) returns the exponential of a chebfun2. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% check for empty chebfun2
if ( isempty( f ) ) 
    return 
end 

op = @(x,y) exp( feval(f, x, y) ); % resample.
f = chebfun2( op, f.domain );      % Call constructor.

end