function [SA,V,zC,I] = revolution(f)
%REVOLUTION   Calculated quantities for a surface of revolution.
%   Given a chebfun F that describes the radius of a surface of revolution,
%   REVOLUTION(F) returns a structure with the surface area, enclosed
%   volume, (nontrivial coordinate of) the center of mass, and the moment
%   of inertia of the solid.
%
%   [SA,V,ZC,I] = REVOLUTION(f) returns the items separately.
%
%   Example:
%
%      >> f = chebfun(@(z) 2-2*z,[0 1]);
%      >> cheb.revolution(f)
%      ans =
%
%             surfaceArea: 14.0496
%                  volume: 4.1888
%               centroidZ: 0.2500
%         momentOfInertia: 5.0265

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if min(f) < -1e-15
    error('Radius of the revolved surface must be nonnegative.')
end

z = chebfun(@(z) z,domain(f));
% Compute the quantities of interest using well known formulas from calculus.
SA = 2*pi*sum(f.*sqrt(1 + abs(diff(f)).^2));
V = pi*sum(f.^2);
zC = pi/V*sum(z.*f.^2);
I = pi/2*sum(f.^4);

if nargout <= 1
    result.surfaceArea = SA;
    result.volume = V;
    result.centroidZ = zC;
    result.momentOfInertia = I;
    SA = result;
end
