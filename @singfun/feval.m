function val = feval(f,x)
%FEVAL Evaluate a singfun at given points.
%   TODO: Documentation

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Evaluate the smooth part.
val = feval(f.smoothPart,x);

% Apply the values from the exponents and scale correctly
val = val.*(x+1).^f.exponents(1);
val = val.*(1-x).^f.exponents(2);

end
