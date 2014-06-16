function pass = test_ivp(pref)
% This tests solves a linear IVP using chebop and checks
% that the computed solution is accurate.

% NOTE: Taken from V4 test chebop_ivp.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

%%

d = [-1, 1];
x = chebfun(@(x) x, d);
A = chebop(@(x, u) diff(u) - u, d);
A.lbc = exp(-1) - 1;
u = mldivide(A, 1-x, pref);
err(1) = norm( u - (exp(x) + x) );

pass = err < tol;

end
