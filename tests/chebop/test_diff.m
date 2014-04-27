function pass = test_diff(pref)
% Checks if the differentiation chebop is equivalent to differentiating
% a chebfun.

% NOTE: Taken from V4 test chebop_diff.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

d = [-3, -1.5];
D = chebop(@(u) diff(u), d);
f = chebfun(@(x) exp(sin(x).^2+2), d);
pass(1) = norm(D*f - diff(f)) < tol;

end
