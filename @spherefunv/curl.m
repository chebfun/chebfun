function G = curl(F) 
%CURL  Numerical surface curl.
%   G = CURL(F) returns the numerical surface curl of the SPHEREFUNV F.
%   This only makes mathematical sense for a tangential vector field. 
%
% See also SPHEREFUNV/DIVERGENCE, SPHEREFUN/GRADIENT, SPHEREFUNV/VORTICITY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check.
if ( isempty(F) )
    G = spherefunv;
    return
end

Fc = F.components;

% The following derivatives are all tangential derivatives.

% First component: d/dy(fz) - d/dz(fy)
gx = diff(Fc{3}, 2) - diff(Fc{2}, 3);

% Second component: d/dz(fx) - d/dx(fz)
gy = diff(Fc{1}, 3) - diff(Fc{3}, 1);

% Third component: d/dx(fy) - d/dy(fx)
gz = diff(Fc{2}, 1) - diff(Fc{1}, 2);

G = spherefunv(gx,gy,gz);

end
