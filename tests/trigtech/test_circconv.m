% Test file for trigtech/circconv

function pass = test_circconv(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Generate a random point to use as test values.
seedRNG(6178);
x = 2 * rand - 1;

%%
% Check operation in the face of empty arguments.

f = testclass.make();
g = testclass.make(@(x) sin(pi*x), [], pref);
pass(1) = (isempty(circconv(f,g)) && isempty(circconv(g,f)));

%%
% Simple checks

f_op = @(x) tanh(5*cos(pi*(x)));
g_op = @(x) ones(size(x));
f = testclass.make(f_op, [], pref);
g = testclass.make(g_op, [], pref);

hfg = circconv(f,g);
hgf = circconv(g,f);
% Answer should be zero since the functions are odd.
approx_fg = feval(hfg,x);
approx_gf = feval(hgf,x);
expected = 0;
err = abs(approx_fg - expected);
tol = 1e2*eps;
pass(2) = err < tol;

f_op = @(x) tanh(cos(pi*(x)));
a = x;
fa_op = @(x) tanh(cos(pi*(x-a)));

% Harder tests
f = testclass.make(f_op, [], pref);
g = circconv(f,f);

% Computed answer at zero:
approx = feval(g,0);
% Expected answer:
expected = sum(f.*f);
err = abs(approx - expected);
tol = 1e2*eps*vscale(g);
pass(3) = err < tol;

% Computed answer at x:
approx = feval(g,x);
% Expected answer:
h = testclass.make(fa_op, [], pref);
expected = sum(f.*h);
err = abs(approx - expected);
tol = 1e2*eps*vscale(g);
pass(4) = err < tol;

end
