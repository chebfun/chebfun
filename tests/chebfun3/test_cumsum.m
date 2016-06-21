function pass = test_cumsum(pref)

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Check cumsum on cube domain
[x, y, z] = cheb.xyz;

pass(1) = norm(cumsum(x) - .5*(x.^2-1)) < tol;
pass(2) = norm(cumsum(x, 1) - .5*(x.^2-1)) < tol;
pass(3) = norm(cumsum(x, 2) - x.*(y+1)) < tol;
pass(4) = norm(cumsum(x, 3) - x.*(z+1)) < tol;

pass(5) = norm(cumsum(y) - y.*(x+1)) < tol;
pass(6) = norm(cumsum(y, 1) - y.*(x+1)) < tol;
pass(7) = norm(cumsum(y, 2) - .5*(y.^2-1)) < tol;
pass(8) = norm(cumsum(y, 3) - y.*(z+1)) < tol;

pass(9) = norm(cumsum(z) - z.*(x+1)) < tol;
pass(10) = norm(cumsum(z, 1) - z.*(x+1)) < tol;
pass(11) = norm(cumsum(z, 2) - z.*(y+1)) < tol;
pass(12) = norm(cumsum(z, 3) - .5*(z.^2-1)) < tol;

% Check cumsum on rectangular box domain
x = chebfun3(@(x,y,z) x, [-1.1  2 -0.2 3 5 6]);
y = chebfun3(@(x,y,z) y, [-1.1  2 -0.2 3 5 6]);
z = chebfun3(@(x,y,z) z, [-1.1  2 -0.2 3 5 6]);

pass(13) = norm(cumsum(x) - .5*(x.^2-(-1.1)^2)) < tol;
pass(14) = norm(cumsum(x, 1) - .5*(x.^2-(-1.1)^2)) < tol;
pass(15) = norm(cumsum(x, 2) - x.*(y+0.2)) < tol;
pass(16) = norm(cumsum(x, 3) - x.*(z-5)) < tol;

pass(17) = norm(cumsum(y) - y.*(x+1.1)) < tol;
pass(18) = norm(cumsum(y, 1) - y.*(x+1.1)) < tol;
pass(19) = norm(cumsum(y, 2) - .5*(y.^2-(-0.2)^2)) < tol;
pass(20) = norm(cumsum(y, 3) - y.*(z-5)) < tol;

pass(21) = norm(cumsum(z) - z.*(x+1.1)) < tol;
pass(22) = norm(cumsum(z, 1) - z.*(x+1.1)) < tol;
pass(23) = norm(cumsum(z, 2) - z.*(y+0.2)) < tol;
pass(24) = norm(cumsum(z, 3) - .5*(z.^2-5^2)) < tol;

% Check that a triple cumsum is a cumsum3
f = sin((x-.1).*(y+.1).*(z+.1));
pass(25) = norm(cumsum(cumsum(cumsum(f), 2), 3) - cumsum3(f)) < tol; 
pass(26) = norm(cumsum(cumsum(cumsum(f, 1), 2), 3) - cumsum3(f)) < tol; 

% Look in one direction and make sure we get the right thing:
f = chebfun(@(x) exp(x));
f3 = chebfun3(@(x,y,z) exp(x));
g = cumsum(f);
g3 = cumsum(f3);
g3x = g3(:,0,.3);
pass(27) = norm(g3x-g) < tol;

end