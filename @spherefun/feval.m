function y = feval(f, x, y, z)
%FEVAL  Evaluate a SPHEREFUN at one or more points.
%   Y = FEVAL(F, LAMBDA, THETA)  evaluates a spherefun F at (LAMBDA, THETA)
%   where LAMBDA and THETA are doubles representing the longitudinal (or 
%   azimuthal) and latitudinal (or elevation) angles.
%
%   Y = FEVAL(F, X, Y, Z)  evaluates a spherefun F at a point (X,Y,Z) in
%   Cartesian cooridnates on the surface of a sphere.
%
% See also SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 3 )      % Spherical coordinates used.
    lambda = x;
    theta = y;
    y = feval@separableApprox(f, lambda, theta);
    
elseif ( nargin == 4 ) % Cartesian coordinates used.
    if ( isnumeric(x) && isnumeric(y) && isnumeric(z) )
        % Convert to spherical coordinates
        [lambda, theta, rad] = cart2sph(x,y,z);
        
        if ( any(rad > (1 + 1e-8)) )
            error('CHEBFUN:SPHEREFUN:FEVAL:pointsNotOnSphere',...
                ['The specified points to evaluate the function do not '...
                'lie sufficiently close to the surface of the '...
                'unit sphere.']);
        end
        % Check latitudinal coordinate system to set the elevation angle
        % appropriately.
        if iscolat(f)
            theta = pi/2 - theta;
        end
        y = feval@separableApprox(f, lambda, theta);
        
    elseif ( strcmp(x, ':') && strcmp(y, ':') && strcmp(z, ':') )
        y = f;
        
    elseif ( strcmp(x, ':') && isnumeric(y) && isnumeric(z) )
        y = [ feval(f, sqrt(1-y.^2-z.^2), y, z) ; 
              feval(f, -sqrt(1-y.^2-z.^2), y, z) ];
        
    elseif ( strcmp(y, ':') && isnumeric(x) && isnumeric(z) )
        y = [ feval(f, x, sqrt(1-x.^2-z.^2), z); 
              feval(f, x, -sqrt(1-x.^2-z.^2), z) ];
        
    elseif ( strcmp(z, ':') && isnumeric(x) && isnumeric(y) )
        y = [ feval(f, x, y, sqrt(1-x.^2-y.^2)); 
              feval(f, x, y, -sqrt(1-x.^2-y.^2)) ];
        
    elseif ( strcmp(x, ':') && strcmp(y, ':') && isnumeric(z) )
        y = chebfun(@(t) feval(f, sqrt(1-z.^2).*cos(t), ...
                              sqrt(1-z.^2).*sin(t), z), [-pi, pi], 'trig');
        
    elseif ( strcmp(x, ':') && strcmp(z, ':') && isnumeric(y) )
        y = chebfun(@(t) feval(f, sqrt(1-y.^2).*cos(t), y, ...
                                 sqrt(1-y.^2).*sin(t)), [-pi, pi], 'trig');
        
    elseif ( strcmp(y, ':') && strcmp(z, ':') && isnumeric(x) )
        y = chebfun(@(t) feval(f, x, sqrt(1-x.^2).*cos(t), ...
                                 sqrt(1-x.^2).*sin(t)), [-pi, pi], 'trig');
        
    else
        error('CHEBFUN:SPHEREFUN:feval:argin',['Unkown input '...
            'feval(%s,%s,%s)',x,y,z]);
    end

end

end