function g = diskfun(f, a, b, varargin)
% DISKFUN is the intersection between a BALLFUN function and a
%   plane
%   DISKFUN(f, lambda, theta, r) is the slice of the BALLFUN function f
%   corresponding to the normal (lambda, theta) at radius r
%   DISKFUN(f, 'x', 'y', r) is the slice of the BALLFUN function f
%   corresponding to the plane X-Y at radius r

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Take a ballfun function and a point on the sphere (r, lambda, theta), 
% Return is in [0, 1], lambda is in [-pi, pi], theta is in [0, pi]
% r the function f on the disk parametrized by the normal 0-(r, lambda,
% theta)
% The point (1-r,0) of the disk corresponds to the point (1-r, lambda, pi/2 + theta)
% The point (-(1-r),0) of the disk corresponds to the point (1-r, pi + lambda, pi/2 -
% theta)

% If a and b are char then compute the corresponding angles lambda and theta
if isletter(a) && isletter(b)
    A = plane2angle(a,b);
    lambda = A(1);
    theta = A(2);
else
    lambda = a;
    theta = b;
end

% Find the radius r
if nargin > 3
    r = varargin{1};
else
    r = 0;
end

[m,n,p] = size(f);

% Resize the number of values
%m = max(m, 20);
%n = max(n, 20);

% If n is odd, make it even
n = n + mod(n,2);

% Points of evaluation on the disk
ListR = chebpts(m, [0 1-r]);
ListLam = pi*trigpts(n);

% Convert the polar points to cartesian points
[Lam,R] = meshgrid(ListLam,ListR);
x = R.*cos(Lam);
y = R.*sin(Lam);

% Changement of coordinates
MrotationZ = [cos(lambda), -sin(lambda), 0;...
              sin(lambda),  cos(lambda), 0;...
              0, 0, 1];
MrotationX = [1, 0, 0;...
              0, cos(theta), -sin(theta);...
              0, sin(theta),  cos(theta)];
Mrotation = MrotationZ*MrotationX;

% Evaluate the function at the points r, lam
G = zeros(m, n);
for i = 1:m
    for j = 1:n
        vec_ij = [x(i,j); y(i,j); r];
        % Rotation of the normal vector by lambda and theta
        vec_0 = Mrotation*vec_ij;
        x_0 = vec_0(1); y_0 = vec_0(2); z_0 = vec_0(3); 
        [lam_0, th_0, r_0] = cart2sph(x_0, y_0, z_0);
        th_0 = pi/2-th_0;
        % Evaluate the function at the new point
        G(i,j) = feval(f, r_0, lam_0, th_0);
    end
end

% Return the diskfun
g = diskfun(real(G));
end

% Take to char a and b and return the angles lambda and theta to get the
% normal to the plane a - b
function A = plane2angle(a,b)
if a == 'x' && b == 'y'
    lambda = 0;
    theta = 0;
elseif a == 'x' && b == 'z'
    lambda = 0;
    theta = pi/2;   
elseif a == 'y' && b == 'z'
    lambda = pi/2;
    theta = pi/2;
else
    error('BALLFUN:SLICE:unknown', ...
                ['Undefined function ''slice'' for input arguments ' ...
                '''%s'' and ''%s''.'], a, b);
end
A = [lambda, theta];
end
