function pass = test_simplify(pref)
% [TODO] Make this test more extensive.

% Get preferences:
if ( nargin < 1 )
    pref = funcheb1.pref;
end

% Set the tolerance:
tol = 10*pref.funcheb1.eps;

%%
% Test on a scalar-valued function:

f = @(x) sin(x);
pref.funcheb1.n = 33;
g = funcheb1(f, 0, pref);
h = simplify(g);
x = funcheb1.chebpts(14);
pass(1) = length(g) == 33 && length(h) == 14 && norm(f(x) - h.values, inf) < tol;

%%
% Test on a vector-valued function:

f = @(x) [ sin(x), cos(x), exp(x) ];
pref.funcheb1.n = 33;
g = funcheb1(f, 0, pref);
h = simplify(g);
x = funcheb1.chebpts(15);
pass(2) = length(g) == 33 && length(h) == 15 && norm(f(x) - h.values, inf) < tol;

end