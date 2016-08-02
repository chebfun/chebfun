function G = grad( f ) 
%GRAD  Numerical gradient of a DISKFUN. 
%   G = GRAD(F) returns the numerical gradient of the
%   DISKFUN F as a DISKFUNV G.
%
%   This is shorthand for the command GRADIENT.
%
% See also GRADIENT, DIV, CURL, VORT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

G = gradient( f );

end

