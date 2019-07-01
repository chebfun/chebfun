function pass = test_points(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

n = 10;
tol = 10*eps;
r = rand(n,1);
dom = [0, 1];

for k = 1:2
    f = chebfun(r, dom, 'chebkind', k);
    pass(k) = norm(f.points - chebpts(n, dom, k)) < tol;
end