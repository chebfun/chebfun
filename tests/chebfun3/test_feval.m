function pass = test_feval( pref ) 
% Test feval

if ( nargin == 0) 
    pref = chebfunpref; 
end

seedRNG(42);

tol = 100*pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(j) = (abs(f(0,0,0)) < tol*vscale(f));  
j = j+1;

pass(j) = (abs(f(pi/6,pi/12,-1)-pi/6) < tol*vscale(f));  
j = j+1;

f = chebfun3(@(x,y,z) y, [-1 2 -pi/2 pi -3 1]); 
pass(j) = (abs(f(0,0,0)) < tol);   
j = j+1;

pass(j) = (abs(f(pi/6,pi/12,-1)-pi/12) < tol*vscale(f)); 
j = j+1;

f = chebfun3(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(j) = (abs(f(0,0,0)) < tol);   
j = j+1;

pass(j) = (abs(f(pi/6,pi/12,-1)+1) < tol*vscale(f)); 
j = j+1;

% some harder tests. 
f = @(x,y,z) cos(x) + sin(x.*y) + sin(z.*x); 
g = chebfun3(f);
pts = 2*rand(3,1) - 1;
pass(j) = (abs(f(pts(1),pts(2),pts(3)) - g(pts(1),pts(2),pts(3)))<tol*vscale(g));
j = j+1;

% Are we evaluating on arrays correctly?
r = rand(10,1); 
s = rand(10,1); 
t = rand(10,1); 
[rr, ss, tt]=meshgrid(r,s,t);
pass(j) = (norm((f(r,s,t) - g(r,s,t))) < tol*vscale(g));
j = j+1;

pass(j) = (max(max(max(abs(f(rr,ss,tt) - g(rr,ss,tt))))) < tol*vscale(g));
j = j+1;

% Does this work off [-1,1]^2
g = chebfun3(f,[-pi/6 pi/2 -pi/12 sqrt(3) -3 1]); % strange domain. 
r = 0.126986816293506; s = 0.632359246225410; t = 0.351283361405006;
% three fixed random number in domain.
pass(j) = (abs(f(r,s,t) - g(r,s,t))<tol*vscale(g));
j = j+1;

% Are we evaluating on arrays correctly
pass(j) = (norm((f(r,s,t) - g(r,s,t)))<tol*vscale(g));
j = j+1;

pass(j) = (max(max(max(abs(f(rr,ss,tt) - g(rr,ss,tt)))))<tol*vscale(g)); 
j = j+1;

% vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
xx = linspace(-1, 1, 100)';
yy = linspace(-1, 1, 100)';
zz = linspace(-1, 1, 100)';
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F - FF) < 100*tol;
j = j+1;

% 'trig' flag + vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff, 'trig');
xx = linspace(-1, 1, 100)';
yy = linspace(-1, 1, 100)';
zz = linspace(-1, 1, 100)';
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F - FF) < 100*tol;
j = j+1;

% random vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
xx = rand(100, 1);
yy = rand(100, 1);
zz = rand(100, 1);
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F - FF) < 100*tol;
j = j+1;

% 'trig' flag + random vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff, 'trig');
xx = rand(100, 1);
yy = rand(100, 1);
zz = rand(100, 1);
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F - FF) < 100*tol;
j = j+1;

% meshgrid
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
[xx, yy, zz] = meshgrid(linspace(-1, 1, 100));
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F(:) - FF(:)) < 100*tol;
j = j+1;

% 'trig' flag + meshgrid
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff, 'trig');
[xx, yy, zz] = meshgrid(linspace(-1, 1, 100));
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F(:) - FF(:)) < 100*tol;
j = j+1;

% ndgrid
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
[xx, yy, zz] = ndgrid(linspace(-1, 1, 100));
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F(:) - FF(:)) < 100*tol;
j = j+1;

% 'trig' flag + ndgrid
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff, 'trig');
[xx, yy, zz] = ndgrid(linspace(-1, 1, 100));
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F(:) - FF(:)) < 100*tol;
j = j+1;

% random tensor inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
xx = rand(10, 20, 30);
yy = rand(10, 20, 30);
zz = rand(10, 20, 30);
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F(:) - FF(:)) < 100*tol;
j = j+1;

% 'trig' flag + random tensor inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff, 'trig');
xx = rand(10, 20, 30);
yy = rand(10, 20, 30);
zz = rand(10, 20, 30);
F = f(xx,yy,zz);
FF = ff(xx,yy,zz);
pass(j) = norm(F(:) - FF(:)) < 100*tol;
j = j+1;
end