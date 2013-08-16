% Test file for singfun/plus.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref.singfun;
end

% Generate a few random points to use as test values.
seedRNG(786);
x = -1 + 2*rand(100, 1);

pass = zeros(1, 12); % Pre-allocate pass vector

%%
% Check feval on a simple SINGFUN
fh = @(x) 1./((1+x).*(1-x));
f = singfun(fh, [-1, -1]);
tol = 1e3*pref.eps;
pass(1) = norm(feval(f,x) - feval(fh,x), inf) < tol;

%%
end