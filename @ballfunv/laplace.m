function f = laplace(v)
% LAPLACE Laplacian of a BALLFUNV
%   LAPLACE(v) is the laplacian of the BALLFUNV v

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = grad(div(v))-curl(curl(v));
end
