function g = diskfun(f)
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

% rotate the sphere then evaluate at theta = pi/2

[m,n,p] = size(f);

% If n is odd, make it even
m = m + 1-mod(m,2);
n = n + mod(n,2);

% Evaluation points in [0,1]
r = chebpts(m);
r = r(ceil(m/2):end);

% Evaluation points in [-pi,pi[
lambda = pi*trigpts(n);

[rr, ll, tt] = ndgrid(r, lambda, pi/2);

G = feval(f,rr,ll,tt);

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
