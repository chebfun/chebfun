function f = cos( f ) 
% COS   Cosine of a chebfun2.
%
% See also COSH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Check for empty
if ( isempty( f ) )
    return
end 

op = @(x,y) cos( feval( f, x, y ) );  % Resample. 
f = chebfun2(op, f.domain);           % Call constructor. 

end