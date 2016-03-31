function f = div(f)
%DIV   Numerical surface divergence of a SPHEREFUNV. 
%   D = DIVERGENCE(F) returns the numerical surface divergence of the
%   SPHEREFUNV. This operations only makes mathematical sense for F that
%   are tanget to the sphere.
%
%  This is shorthand for the command DIVERGENCE. 
% 
% See also SPHEREFUNV/DIVERGENCE, SPHEREFUN/GRAD, SPHEREFUNV/CURL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


f = divergence(f); 

end
