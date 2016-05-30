function pass = test_constructor(pref) 
% This tests the chebfun3 constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;

% Can we make a chebfun3: 
ff = @(x,y,z) cos(x+z) + sin(x.*y.*z);  % simple function
fstr = 'cos(x+z) + sin(x.*y.*z)';       % string version

% Construct with the default domain:
f1 = chebfun3(ff); 
f2 = chebfun3(ff, [-1 1 -1 1 -1 1]);
pass(1) = norm(f1 - f2) < tol; 

% Construct from a string:
f3 = chebfun3(fstr);
f4 = chebfun3(fstr, [-1 1 -1 1 -1 1]);
pass(2) = norm(f1 - f3) < tol;
pass(3) = norm(f3 - f4) < tol;

f5 = chebfun3(@(x,y,z) f4(x,y,z));
pass(4) = norm(f3 - f5) < tol;

f5 = chebfun3(@(x,y,z) f4(x,y,z), 'fiberDim', 1);
pass(5) = norm(f3 - f5) < tol;

f5 = chebfun3(@(x,y,z) f4(x,y,z), 'fiberDim', 2);
pass(6) = norm(f3 - f5) < tol;

f5 = chebfun3(@(x,y,z) f4(x,y,z), 'fiberDim', 3);
pass(7) = norm(f3 - f5) < tol;

% Check accuracy:
g = @(x,y,z) cos(x).*y.*z + x.*z.*sin(y); 
f = chebfun3(g);
pass(8) = abs(feval(f, .1, pi/6, -0.3) - feval(g, .1, pi/6, -0.3) ) < tol;

% Check that multilinear rank is computed correct.
g = @(x,y,z) x.*y.*z + y.^2.*z;
f = chebfun3(g);
[r1, r2, r3] = rank(f); %The x-rank should be 2, the y-rank should be 2 and
% the z-rank should be 1.
pass(9) = r1 == 2;
pass(10) = r2 == 2;
pass(11) = r3 == 1;

f = @(x,y,z) 1./(1+x.^2.*y.^2.*z.^2);
ffch = chebfun3(@(x,y,z) f(x,y,z), [-2 2 -2 2 -2 2]);
xx = linspace(-2,2); 
[XX,YY, ZZ] = meshgrid(xx);
pass(12) = max(max(max(abs(f(XX,YY,ZZ) - ffch(XX,YY,ZZ))))) < 1e5*tol;


f = chebfun3(@(x,y,z) 1, 'vectorize');
g = chebfun3(1);
pass(13) = norm(f - g) < tol;
f = chebfun3(@(x,y,z) 1, 'vectorize');
g = chebfun3(1);
pass(14) = norm(f - g) < tol;

% Test length of a chebfun3 in different directions:
f = chebfun3(@(x,y,z) sin(80*x+y+z));
pass(15) = length(f.rows(:,1)) < 50;
pass(16) = length(f.tubes(:,1)) < 50;

% Test whether there is a huge overesimation in the length of a chebfun3:
f = chebfun(@(x) 1./(1+25*x.^2));
f2 = chebfun3(@(x,y,z) 1./(1+25*x.^2));
pass(17) = length(f2.rows) < length(f)+20;

end