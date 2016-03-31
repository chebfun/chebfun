function G = grad(f) 
%GRAD   Numerical surface gradient of a SPHEREFUN. 
%   G = GRAD(F) returns the numerical surface gradient of the SPHEREFUN F as a
%   SPHEREFUNV G.
%
%   This is shorthand for the command GRADIENT.
%
% See also SPHEREFUN/GRADIENT, SPHEREFUNV/DIV, SPHEREFUNV/CURL, SPHEREFUNV/VORT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

G = gradient(f);

end

