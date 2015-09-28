function y = feval( f, c1, c2, c3)
%FEVAL  Evaluate a SPHEREFUN at one or more points.
%   Y = FEVAL( F, LAMBDA, THETA )  evaluates a spherefun F at (LAMBDA,
%   THETA) where LAMBDA and THETA are doubles representing the longitudinal
%   (or azimuthal) and latitudinal (or elevation) angles.
%
%   Y = FEVAL( F, X, Y, Z )  evaluates a spherefun F at a point (X,Y,Z) in
%   Cartesian cooridnates on the surface of a sphere.  
%
% See also SUBSREF.

if nargin == 3      % Spherical coordinates used.
    lambda = c1;
    theta = c2;
elseif nargin == 4  % Cartesian coordinates used.
    if ~isnumeric(c1) || ~isnumeric(c1) || ~isnumeric(c1)
        if strcmpi(c1, ':') || strcmpi(c2, ':') || strcmpi(c3, ':') 
            error('SPHEREFUN:feval:colon','Colon operator not allowed when using Cartesian coordinates');
        else
            error('SPHEREFUN:feval:unknown','Unkown input feval(%s,%s,%s)',c1,c2,c3);
        end
    end
    % Convert to spherical coordinates
    [lambda,theta] = cart2sph(c1,c2,c3);

    % Check latitudinal coordinate system to set the elevation angle
    % appropriately.
    if iscolat( f )
        theta = pi/2 - theta;
    end
end

y = feval@separableApprox( f, lambda, theta );

if size(lambda,1) == 1 && size(theta,1) == 1
    y = y.';
end

end 