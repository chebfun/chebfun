% Test file for funcheb2/max.

function pass = test_max(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

%%
% Spot-check the extrema for a few functions.

pass(1) = test_spotcheck_max(@(x) ((x-0.2).^3 -(x-0.2) + 1).*sec(x-0.2), ...
    1.884217141925336, pref);
pass(2) = test_spotcheck_max(@(x) sin(10*x), 1, pref);
pass(3) = test_spotcheck_max(@airy, airy(-1), pref);
pass(4) = test_spotcheck_max(@(x) -1./(1 + x.^2), -0.5, pref);
pass(5) = test_spotcheck_max(@(x) (x - 0.25).^3.*cosh(x), 0.75^3*cosh(1), pref);

%%
% Check operation for vectorized inputs.

fun_op = @(x) [sin(10*x) airy(x) (x - 0.25).^3.*cosh(x)];
f = funcheb2(fun_op, [], pref);
[y, x] = max(f);
exact_max = [1 airy(-1) 0.75^3*cosh(1)];
fx = [sin(10*x(1)) airy(x(2)) (x(3) - 0.25).^3.*cosh(x(3))];
pass(6) = (all(abs(y - exact_max) < 10*f.epslevel) && ...
           all(abs(fx - exact_max) < 10*f.epslevel));

% TODO:  Test for complex-valued funcheb2 objects?

end

% Spot-check the results for a given function.
function result = test_spotcheck_max(fun_op, exact_max, pref)

f = funcheb2(fun_op, [], pref);
[y, x] = max(f);
fx = fun_op(x);
result = ((abs(y - exact_max) < 10*f.epslevel) && ... 
          (abs(fx - exact_max) < 10*f.epslevel));

end
