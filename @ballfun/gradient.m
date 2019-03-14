function v = gradient(f)
%GRADIENT Gradient of a BALLFUN in cartesian coordinates.
%   GRADIENT(F) is the gradient of the BALLFUN F expressed in
%   cartesian coordinates.
%
% See also DIV, CURL

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Vx = diff(f, 1);
Vy = diff(f, 2);
Vz = diff(f, 3);
v = ballfunv(Vx, Vy, Vz);
end