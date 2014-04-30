function pass = test_chebop_diff
% Checks if the differentiation chebop is equivalent to differentiating
% a chebfun.

p = chebpref;

d = [-3,-1.5];
D = chebop(@(u) diff(u), d);
f = chebfun(@(x) exp(sin(x).^2+2),d);
pass(1) = norm(D*f - diff(f)) < p.eps;

D2 = 2*D;
f = chebfun(@(x) exp(sin(x).^2+2),d);
pass(2) = norm(D2*f - 2*diff(f)) < p.eps;
