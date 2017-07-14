function pass = test_ultra2ultra(pref)
% Test for ultra2ultra

% TO test, we construct a chebfun and obtain it's coefficients using
% chebfun.ultracoeffs for two values of lamba. We then convert between these
% expansions using ultra2ultra and check the difference is within tolerance. 
% (This test therefore assumes that chebfun.ultracoeffs is behaving correctly.)

if ( nargin == 0 )
    pref = chebfunpref;
end

tol = 100*eps;

%% Real-valued function

f = chebfun(@exp, pref);
lam1 = .6;
lam2 = .7;
c1 = ultracoeffs(f, lam1);
c2 = ultracoeffs(f, lam2);
pass(1) = norm(ultra2ultra(c1, lam1, lam2) - c2) < tol;
pass(2) = norm(ultra2ultra(c2, lam2, lam1) - c1) < tol;

%% Complex-valued function

f = chebfun(@(x) exp(1i*x), [0 3], pref);
lam1 = .6;
lam2 = .7;
c1 = ultracoeffs(f, lam1);
c2 = ultracoeffs(f, lam2);
pass(3) = norm(ultra2ultra(c1, lam1, lam2) - c2) < tol;
pass(4) = norm(ultra2ultra(c2, lam2, lam1) - c1) < tol;

end