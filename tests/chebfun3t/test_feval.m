function pass = test_feval(pref)
% Test CHEBFUN3/FEVAL

if ( nargin == 0)
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

seedRNG(42);

f = chebfun3t(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(1) = abs(f(0,0,0)) < tol*f.vscale;

pass(2) = abs(f(pi/6,pi/12,-1)-pi/6) < tol*f.vscale;  

f = chebfun3t(@(x,y,z) y, [-1 2 -pi/2 pi -3 1]); 
pass(3) = abs(f(0,0,0)) < tol;

pass(4) = abs(f(pi/6,pi/12,-1)-pi/12) < tol*f.vscale;

f = chebfun3t(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(5) = abs(f(0,0,0)) < tol;

pass(6) = abs(f(pi/6,pi/12,-1)+1) < tol*f.vscale;

% some harder tests. 
f = @(x,y,z) cos(x) + sin(x.*y) + sin(z.*x);
g = chebfun3t(f);
pts = 2*rand(3,1) - 1;
pass(7) = abs(f(pts(1),pts(2),pts(3)) - g(pts(1),pts(2),pts(3))) < ...
    tol*g.vscale;

% Are we evaluating on arrays correctly?
r = rand(10,1); 
s = rand(10,1); 
t = rand(10,1); 
[rr, ss, tt]=meshgrid(r,s,t);
pass(8) = max(abs((f(r,s,t) - g(r,s,t)))) < tol*g.vscale;

pass(9) = max(max(max(abs(f(rr,ss,tt) - g(rr,ss,tt))))) < tol*g.vscale;

% Does this work off [-1,1]^2
g = chebfun3t(f,[-pi/6 pi/2 -pi/12 sqrt(3) -3 1]); % strange domain. 
r = 0.126986816293506; s = 0.632359246225410; t = 0.351283361405006;
% three fixed random number in domain.
pass(10) = abs(f(r,s,t) - g(r,s,t))<tol*g.vscale;

% Are we evaluating on arrays correctly
pass(11) = abs(f(r,s,t) - g(r,s,t)) < tol*g.vscale;

pass(12) = max(max(max(abs(f(rr,ss,tt) - g(rr,ss,tt))))) < tol*g.vscale;

% vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
xx = linspace(-1, 1, 100)';
yy = linspace(-1, 1, 100)';
zz = linspace(-1, 1, 100)';
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(13) = norm(F(:) - FF(:)) < 100*tol;

% random vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
xx = rand(100, 1);
yy = rand(100, 1);
zz = rand(100, 1);
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(14) = norm(F - FF) < 100*tol;

% meshgrid
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
[xx, yy, zz] = meshgrid(linspace(-1, 1, 100));
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(15) = norm(F(:) - FF(:)) < 100*tol;

% ndgrid
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
[xx, yy, zz] = ndgrid(linspace(-1, 1, 100));
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(16) = norm(F(:) - FF(:)) < 100*tol;

% random tensor inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
xx = rand(10, 20, 30);
yy = rand(10, 20, 30);
zz = rand(10, 20, 30);
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(17) = norm(F(:) - FF(:)) < 100*tol;

end