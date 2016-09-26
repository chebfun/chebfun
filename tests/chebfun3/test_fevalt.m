function pass = test_fevalt(pref)
% Test fevalt

if ( nargin == 0 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb3Prefs.chebfun3eps;

ff = @(x,y,z) cos(x) + sin(x.*y) + sin(z.*x); 
f = chebfun3(ff);
seedRNG(42);
x = rand(10,1);
y = rand(10,1);
z = rand(10,1);
F = fevalt(f,x,y,z);

[xx, yy, zz] = ndgrid(x,y,z);
Fexact = ff(xx,yy,zz);
pass(1) = norm(Fexact(:) - F(:)) < tol*vscale(f);

dom = [-1 1 -pi/3 pi/3 -2 0];
f = chebfun3(ff, dom);
x = chebpts(20,dom(1:2));
y = chebpts(20,dom(3:4));
z = chebpts(20,dom(5:6));
F = fevalt(f,x,y,z);

[xx, yy, zz] = ndgrid(x,y,z);
Fexact = ff(xx,yy,zz);
pass(2) = norm(Fexact(:) - F(:)) < tol*vscale(f);

% 'trig' flag + vector inputs
ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff, 'trig');
x = linspace(-1, 1, 20)';
y = linspace(-1, 1, 20)';
z = linspace(-1, 1, 20)';
F = fevalt(f,x,y,z);

[xx,yy,zz] = ndgrid(x,y,z);
Fexact = ff(xx,yy,zz);
pass(3) = norm(F(:) - Fexact(:)) < tol*vscale(f);

end