function pass = test_constructor(pref) 
% This tests the chebfun3t constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;

% Can we make a chebfun3t?
ff = @(x,y,z) cos(x+z) + sin(x.*y.*z);  % simple function
fstr = 'cos(x+z) + sin(x.*y.*z)';       % string version

% Construct with the default domain:
f1 = chebfun3t(ff); 
f2 = chebfun3t(ff, [-1 1 -1 1 -1 1]);

[xx, yy, zz] = meshgrid(linspace(-1,1));
pass(1) = max(max(max(abs(f1(xx,yy,zz) - ff(xx,yy,zz))))) < tol;
pass(2) = max(max(max(abs(f2(xx,yy,zz) - ff(xx,yy,zz))))) < tol;
pass(3) = max(max(max(abs(f1(xx,yy,zz) - f2(xx,yy,zz))))) < tol;

% Construct from a string:
f3 = chebfun3t(fstr);
f4 = chebfun3t(fstr, [-1 1 -1 1 -1 1]);
pass(4) = max(max(max(abs(f1(xx,yy,zz) - f3(xx,yy,zz))))) < tol;
pass(5) = max(max(max(abs(f3(xx,yy,zz) - f4(xx,yy,zz))))) < tol;

% Check accuracy:
g = @(x,y,z) cos(x).*y.*z + x.*z.*sin(y); 
f = chebfun3t(g);
pass(6) = abs(feval(f, .1, pi/6, -0.3) - feval(g, .1, pi/6, -0.3) ) < tol;

f = @(x,y,z) 1./(1+x.^2.*y.^2.*z.^2);
ffch = chebfun3t(@(x,y,z) f(x,y,z), [-2 2 -2 2 -2 2]);
[xx,yy, zz] = meshgrid(linspace(-2,2));
pass(7) = max(max(max(abs(f(xx,yy,zz) - ffch(xx,yy,zz))))) < tol;

% Test length of a chebfun3t in different directions:
f = chebfun3t(@(x,y,z) sin(80*x+y+z));
pass(8) = size(f.coeffs, 2) < 50;
pass(9) = size(f.coeffs, 3) < 50;

% Make sure there is no huge overesimation in the length of chebfun3t objects:
f = chebfun(@(x) 1./(1+25*x.^2));
f2 = chebfun3t(@(x,y,z) 1./(1+25*x.^2));
pass(10) = size(f2.coeffs, 2) < length(f) + 20;

end