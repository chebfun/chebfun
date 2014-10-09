% Test file for trigtech/max.m

function pass = test_max(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
% Spot-check the extrema for a few functions.
pass(1) = test_spotcheck_max(testclass, @(x) exp(-cos(2*pi*x)), exp(1), pref);
pass(2) = test_spotcheck_max(testclass, @(x) sin(10*pi*x), 1, pref);
pass(3) = test_spotcheck_max(testclass, @(x) exp(sin(pi*x).^100), exp(1), pref);
pass(4) = test_spotcheck_max(testclass, @(x) exp(-sin(pi*x).^100), 1, pref);
% Approx to sign function
pass(5) = test_spotcheck_max(testclass, @(x) 4/pi*(sin(pi*x) + ...
    1/3*sin(3*pi*x) + 1/5*sin(5*pi*x) + 1/7*sin(7*pi*x) + 1/9*sin(9*pi*x)), 1.182328208857607, pref);

%%
% Check operation for array-valued inputs.
fun_op = @(x) [exp(-cos(2*pi*x)) sin(10*pi*x) exp(-sin(pi*(x-0.32)).^100)];
f = testclass.make(fun_op, [], pref);
[y, x] = max(f);
exact_max = [exp(1) 1 1];
fx = [exp(-cos(2*pi*x(1))) sin(10*pi*x(2)) exp(-sin(pi*(x(3)-0.32)).^100)];
pass(6) = (all(abs(y - exact_max) < 10*f.epslevel) && ...
           all(abs(fx - exact_max) < 10*f.epslevel));

%%
% Test for complex-valued TRIGTECH objects.
pass(7) = test_spotcheck_max(testclass, ...
    @(x) cos(pi*x) + exp(1i*pi*x), ...
    2, pref);

end

% Spot-check the results for a given function.
function result = test_spotcheck_max(testclass, fun_op, exact_max, pref)

f = testclass.make(fun_op,[], pref);
[y, x] = max(f);
fx = fun_op(x);
result = (all(abs(y - exact_max) < 10*f.vscale.*f.epslevel) && ...
          all(abs(fx - exact_max) < 10*f.vscale.*f.epslevel));

end
