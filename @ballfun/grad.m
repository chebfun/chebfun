function v = grad(f)
%GRAD Gradient of a BALLFUN in cartesian coordinates.
%   GRAD(F) is the gradient of the BALLFUN F expressed in
%   cartesian coordinates.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Vx = diff(f, 1, "cart");
Vy = diff(f, 2, "cart");
Vz = diff(f, 3, "cart");
v = ballfunv(Vx, Vy, Vz);
end