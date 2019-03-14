function f = div(v)
%DIV   Numerical divergence of a BALLFUNV. 
%   D = DIV( F ) returns the numerical divergence of the
%   BALLFUNV. 
%
%  This is shorthand for the command DIVERGENCE. 
% 
% See also DIVERGENCE, GRAD, CURL.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = divergence(v);
end