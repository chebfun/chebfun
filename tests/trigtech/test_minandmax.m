% Test file for trigtech/minandmax.m

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();
%%
% Spot-check the extrema for a few functions.

pass(1) = test_spotcheck_minmax(testclass, @(x) exp(-cos(2*pi*x)), ...
    exp(-1), exp(1), pref);
pass(2) = test_spotcheck_minmax(testclass, @(x) sin(10*pi*x), -1, 1, pref);
pass(3) = test_spotcheck_minmax(testclass, @(x) exp(sin(pi*x).^100), ...
    1, exp(1), pref);
pass(4) = test_spotcheck_minmax(testclass, @(x) exp(-sin(pi*x).^100), ...
    exp(-1), 1, pref);
pass(5) = test_spotcheck_minmax(testclass, @(x) 4/pi*(sin(pi*x) + ...
    1/3*sin(3*pi*x) + 1/5*sin(5*pi*x) + 1/7*sin(7*pi*x) + ...
    1/9*sin(9*pi*x)), -1.182328208857607, 1.182328208857607, pref);

%%
% Check operation for array-valued inputs.

fun_op = @(x) [exp(-cos(2*pi*x)) sin(10*pi*x) exp(-sin(pi*(x-0.32)).^100)];
f = testclass.make(fun_op, [], pref);
[y, x] = minandmax(f);
y_exact = [exp(-1) -1 exp(-1);
           exp(1)   1 1];

pass(6) = all(abs(y(:) - y_exact(:)) < 100*max(vscale(f)*eps));

% Check that the points x are indeed extreme points of the function 
% operator.
for k = 1:1:size(f.coeffs, 2)
    fx = fun_op(x(:, k));
    if ( max(abs(fx(:, k) - y_exact(:, k))) > 10*max(vscale(f)*eps) )
        pass(6) = 0;
        break;
    end
end

% Test complex-array-valued TRIGTECH objects.
f = testclass.make(@(x) [exp(sin(2*pi*x)), 1i*cos(20*pi*x)]);
[vals, pos] = minandmax(f);
f1 = testclass.make(@(x) exp(sin(2*pi*x)));
[vals1, pos1] = minandmax(f1);
f2 = testclass.make(@(x) 1i*cos(20*pi*x));
[vals2, pos2] = minandmax(f2);
pass(7) = norm(abs(vals) - abs([vals1 vals2]), inf) < ...
    1e2*max(vscale(f)*eps);
end

% Spot-check the results for a given function.
function result = test_spotcheck_minmax(testclass, fun_op, exact_min, ...
    exact_max, pref)

f = testclass.make(fun_op, [], pref);
[y, x] = minandmax(f);
y_exact = [exact_min ; exact_max];
fx = fun_op(x);
result = ((max(abs(y - y_exact)) < 100*vscale(f)*eps) && ... 
          (max(abs(fx - y_exact)) < 10*vscale(f)*eps));
    
end
