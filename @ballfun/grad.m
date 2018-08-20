function v = grad(f)
% GRAD Gradient of a BALLFUN function in cartesian system
%   GRAD(f) is the gradient of the BALLFUN function f expressed in the
%   cartesian system

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Vx = diff(f,1,"cart");
Vy = diff(f,2,"cart");
Vz = diff(f,3,"cart");
v = ballfunv(Vx,Vy,Vz);
end