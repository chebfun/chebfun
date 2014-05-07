function pass = test_cumsum(pref)
% Check if the indefinite integral chebop works.

% NOTE: Taken from V4 test chebop_cumsum

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

d = [4, 5.6];
Q = chebop(@(u) cumsum(u), d);
f = chebfun(@(x) exp(sin(x).^2+2), d);
pass(1) = norm(Q*f - cumsum(f)) < tol;

[L, res, isLinear] = linearize(Q);
pass(2) = all(isLinear);

end
