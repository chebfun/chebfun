function y = feval(f,x)
%FEVAL Evaluate a singfun at given points.
%   TODO: Documentation

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Evaluate the smooth part.
y = feval(f.smoothPart,x);

% Apply the values from the exponents and scale correctly
y = y.*(x+1).^f.exponents(1);
y = y.*(1-x).^f.exponents(2);

end
