function f = randnfunball(lambda)
%RANDNFUNBALL   Smooth random function on the unit ball
%   F = RANDNFUNBALL(LAMBDA) returns a BALLFUN of maximum
%   frequency about 2pi/LAMBDA and standard normal distribution N(0,1)
%   at each point.  F is obtained by restricting output of RANDNFUN3
%   to the unit ball.
%
%   RANDNFUNBALL() uses the default value LAMBDA = 1.
%
% See also RANDNFUN3, RANDNFUNSPHERE, RANDNFUNDISK.

% Copyright 2018 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 0
    lambda = 1;
end

% Create a randnfun3 function
fcart = randnfun3( lambda );

% Create a tensor product grid in (r,lam,th)
sz = length(fcart);

m = sz + 1- mod(sz,2);
n = 2*sz;
p = 2*sz;

% Sample on a [sz,sz,sz] grid:
r = chebpts( m );
lam  = pi*trigpts( n );
th  = [pi*trigpts( p ); pi];

[rr, ll, tt] = ndgrid(r(ceil(m/2):end), lam, th(p/2+1:end));

% Sample the random function over the square and create the ballfun:
xx = rr.*sin(tt).*cos(ll);
yy = rr.*sin(tt).*sin(ll);
zz = rr.*cos(tt);
F = feval(fcart, xx, yy, zz);

% Simplify and return the output
f = ballfun(F);
end