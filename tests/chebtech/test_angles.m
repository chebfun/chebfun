function pass = test_angles(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

n = 10;
tol = 1e-15;

% chebtech1:
x = chebtech1.chebpts(n);
t = chebtech1.angles(n);
pass(1) = norm(cos(t) - x, inf) < tol;

% chebtech2:
x = chebtech2.chebpts(n);
t = chebtech2.angles(n);
pass(2) = norm(cos(t) - x, inf) < tol;

end
