function f = exp(f) 
%EXP  Exponential of a SEPARABLEAPPROX
%   EXP(F) returns the exponential of a SEPARABLEAPPROX. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) ) 
    return 
end 

f = compose( f, @exp ); 

end
