function pass = test_BMCIII( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Example 1
f = @(r,lam,th)r.*cos(lam).*sin(th);
m = 3; n = 4; p = 4;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(1) = norm(F(:)-G(:)) < tol;

% Example 2
f = @(r,lam,th)r.*cos(lam).*sin(th);
m = 5; n = 6; p = 8;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(2) = norm(F(:)-G(:)) < tol;

% Example 3
f = @(r,lam,th)r.*cos(lam).*sin(th);
m = 3; n = 4; p = 6;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(3) = norm(F(:)-G(:)) < tol;

% Example 4
f = @(r,lam,th)r.*sin(lam).*sin(th);
m = 3; n = 4; p = 4;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(4) = norm(F(:)-G(:)) < tol;

% Example 5
f = @(r,lam,th)r.*sin(lam).*sin(th);
m = 5; n = 6; p = 8;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(5) = norm(F(:)-G(:)) < tol;

% Example 6
f = @(r,lam,th)r.*sin(lam).*sin(th);
m = 3; n = 4; p = 6;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(6) = norm(F(:)-G(:)) < tol;

% Example 7
f = @(r,lam,th)r.*cos(th)+0*lam;
m = 3; n = 4; p = 4;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(7) = norm(F(:)-G(:)) < tol;

% Example 8
f = @(r,lam,th)r.*cos(th)+0*lam;
m = 5; n = 6; p = 8;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(8) = norm(F(:)-G(:)) < tol;

% Example 9
f = @(r,lam,th)r.*cos(th)+0*lam;
m = 3; n = 4; p = 6;
F = ballfun.BMCIII(f,m,n,p);
r = reshape(chebpts(m),m,1,1);
lam = reshape(pi*trigpts(n),1,n,1);
th = reshape(pi*trigpts(p),1,1,p);
G = f(r,lam,th);
pass(9) = norm(F(:)-G(:)) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end