function pass = test_chebpolyval(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% NOTE: Since CHEBPOLYVAL() is basicaly a wrapper to CHEBTECH/FEVAL(), this is
% just a simple test.

% Init:
seedRNG(42)
tol = 10*pref.eps;

% Use CHEBPOLYVAL():
n = 10;
c = rand(10, 2);
x = rand(3);
fx1 = chebpolyval(c, x);

% Use CHEBPOLY():
T = chebpoly(0:9);
fx2 = feval(T*flipud(c), x);

% Error:
err = norm(fx1 - fx2, inf);
pass = err < n*tol;

end