function pass = test_chebpolyval(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% NOTE: Since CHEBPOLYVAL() is basically a wrapper to CHEBTECH/FEVAL(), this is
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
Tc = T*flipud(c);
fx2 = feval(Tc, x);

% Error:
err = norm(fx1 - fx2, inf);
pass(1) = err < n*tol;

%%

x = chebfun('x');
f = chebpolyval(c, x);
pass(2) = norm(f - Tc, inf) < n*tol;

end
