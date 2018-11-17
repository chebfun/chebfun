function f = laplacian(v)
% LAPLACIAN Laplacian of a BALLFUNV.
%   LAPLACIAN(V) is the laplacian of the BALLFUNV V.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = grad(div(v))-curl(curl(v));
end
