function pass = test_jaccoeffs(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(42)

% Smooth domain:
n = 10;
a = .1;
b = -.3;
J = jacpoly(0:n, a, b);
c = rand(n+1,1);
f = J*c;
err = norm(c - jaccoeffs(f, n+1, a, b).', inf);
tol = 10*epslevel(f);
pass(1) = err < tol;

% Piecewise domain:
n = 11;
a = 1.1;
b = -.9;
J = jacpoly(0:n, a, b, [-1 0 1]);
c = rand(n+1,1);
f = J*c;
err = norm(c - jaccoeffs(f, n+1, a, b).', inf);
tol = 10*vscale(f)*epslevel(f);
pass(2) = err < tol;

end
