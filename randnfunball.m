function f = randnfunball(lambda,S)
%RANDNFUNBALL   Smooth random function on the unit ball
%   F = RANDNFUNBALL(LAMBDA) returns a BALLFUN of maximum
%   frequency about 2pi/LAMBDA and standard normal distribution N(0,1)
%   at each point.  F is obtained by restricting output of RANDNFUN3
%   to the unit ball.
%
%   RANDNFUNBALL() uses the default value LAMBDA = 1.

% Extract the discretization
m = S(1); n = S(2); p = S(3);

% Create a randnfun3 function
fcart = randnfun3( lambda );

% Create a tensor product grid in (r,lam,th)
sz = length(fcart);

if ( sz > min([m, n, p]) )
    
    % Sample on an [m,n,p] grid:
    r = chebpts( m );
    lam  = pi*trigpts( n );
    th  = pi*trigpts( p );
    [rr, ll, tt] = ndgrid(r, lam, th);
    xx = rr.*sin(tt).*cos(ll);
    yy = rr.*sin(tt).*sin(ll);
    zz = rr.*cos(tt);
    F = feval(fcart, xx, yy, zz);
    f = ballfun(F, 'vals');
    
else
    % Sample on a [sz,sz,sz] grid:
    r = chebpts( sz );
    lam  = pi*trigpts( sz );
    th  = pi*trigpts( sz );
    [rr, ll, tt] = ndgrid(r, lam, th);
    
    % Sample the random function over the square and create the ballfun:
    xx = rr.*sin(tt).*cos(ll);
    yy = rr.*sin(tt).*sin(ll);
    zz = rr.*cos(tt);
    F = feval(fcart, xx, yy, zz);
    ftemp = ballfun(F, 'vals');
    
    % Pad coefficient to make f live at the right discretization size:
    X = zeros(m, n, p);
    cfs = ftemp.coeffs;
    %X(1:sz, ceil((n - sz)/2)+(1-mod(n-sz,2)):n-ceil((n - sz)/2),...
       % ceil((p - sz)/2)+(1-mod(p-sz,2)):p-ceil((p - sz)/2) ) = cfs;
    X(1:sz,floor(n/2)+1-floor(sz/2):floor(n/2)+sz-floor(sz/2),floor(p/2)+1-floor(sz/2):floor(p/2)+sz-floor(sz/2)) = cfs;
    f = ballfun( X );
end

end