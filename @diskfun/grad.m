function G = grad( f ) 
%GRAD  Numerical gradient of a DISKFUN. 
%   G = GRAD(F) returns the numericalgradient of the
%   DISKFUN F as a DISKFUNV G.
%
%   This is shorthand for the command GRADIENT.
%
% See also GRADIENT, DIV, CURL, VORT

G = gradient( f );

end

