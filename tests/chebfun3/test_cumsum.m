function pass = test_cumsum( pref )

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 100*pref.cheb3Prefs.chebfun3eps;
j = 1; 

% Check cumsum on cube domain
[x, y, z] = cheb.xyz;

pass(j) = norm(cumsum(x) - .5*(x.^2-1)) < tol;    j = j + 1;
pass(j) = norm(cumsum(x, 1) - .5*(x.^2-1)) < tol; j = j + 1;
pass(j) = norm(cumsum(x, 2) - x.*(y+1)) < tol;    j = j + 1;
pass(j) = norm(cumsum(x, 3) - x.*(z+1)) < tol;    j = j + 1;

pass(j) = norm(cumsum(y) - y.*(x+1)) < tol;       j = j + 1;
pass(j) = norm(cumsum(y, 1) - y.*(x+1)) < tol;    j = j + 1;
pass(j) = norm(cumsum(y, 2) - .5*(y.^2-1)) < tol; j = j + 1;
pass(j) = norm(cumsum(y, 3) - y.*(z+1)) < tol;    j = j + 1;

pass(j) = norm(cumsum(z) - z.*(x+1)) < tol;       j = j + 1;
pass(j) = norm(cumsum(z, 1) - z.*(x+1)) < tol;    j = j + 1;
pass(j) = norm(cumsum(z, 2) - z.*(y+1)) < tol;    j = j + 1;
pass(j) = norm(cumsum(z, 3) - .5*(z.^2-1)) < tol; j = j + 1;

% Check cumsum on rectangular box domain
x = chebfun3(@(x,y,z) x, [-1.1  2 -0.2 3 5 6]);
y = chebfun3(@(x,y,z) y, [-1.1  2 -0.2 3 5 6]);
z = chebfun3(@(x,y,z) z, [-1.1  2 -0.2 3 5 6]);

pass(j) = norm(cumsum(x) - .5*(x.^2-(-1.1)^2)) < tol;    j = j + 1;
pass(j) = norm(cumsum(x, 1) - .5*(x.^2-(-1.1)^2)) < tol; j = j + 1;
pass(j) = norm(cumsum(x, 2) - x.*(y+0.2)) < tol;         j = j + 1;
pass(j) = norm(cumsum(x, 3) - x.*(z-5)) < tol;           j = j + 1;

pass(j) = norm(cumsum(y) - y.*(x+1.1)) < tol;            j = j + 1;
pass(j) = norm(cumsum(y, 1) - y.*(x+1.1)) < tol;         j = j + 1;
pass(j) = norm(cumsum(y, 2) - .5*(y.^2-(-0.2)^2)) < tol; j = j + 1;
pass(j) = norm(cumsum(y, 3) - y.*(z-5)) < tol;           j = j + 1;

pass(j) = norm(cumsum(z) - z.*(x+1.1)) < tol;            j = j + 1;
pass(j) = norm(cumsum(z, 1) - z.*(x+1.1)) < tol;         j = j + 1;
pass(j) = norm(cumsum(z, 2) - z.*(y+0.2)) < tol;         j = j + 1;
pass(j) = norm(cumsum(z, 3) - .5*(z.^2-5^2)) < tol;      j = j + 1;

% Check that a triple cumsum is a cumsum3
f = sin((x-.1).*(y+.1).*(z+.1));
pass(j) = norm(cumsum(cumsum(cumsum(f), 2), 3) - cumsum3(f)) < tol; 
j = j + 1;
pass(j) = norm(cumsum(cumsum(cumsum(f, 1), 2), 3) - cumsum3(f)) < tol; 

end