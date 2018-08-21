function f = divergence(v)
%DIVERGENCE  Divergence of a BALLFUNV in cartesian coordinates.
%   DIVERGENCE(f) is the divergence of the BALLFUNV v expressed in
%   cartesian coordinates.
%
% See also CURL. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[Vx,Vy,Vz] = v.comp{:};
f = diff(Vx,1,"cart") + diff(Vy,2,"cart") + diff(Vz,3,"cart");

end