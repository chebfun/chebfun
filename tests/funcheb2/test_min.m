% Test file for funcheb2/min.

function pass = test_min(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

%%
% Spot-check the extrema for a few functions.

pass(1) = test_spotcheck_min(@(x) -((x-0.2).^3 -(x-0.2) + 1).*sec(x-0.2), ...
    -1.884217141925336, pref);
pass(2) = test_spotcheck_min(@(x) -sin(10*x), -1, pref);
pass(3) = test_spotcheck_min(@(x) -airy(x), -airy(-1), pref);
pass(4) = test_spotcheck_min(@(x) 1./(1 + x.^2), 0.5, pref);
pass(5) = test_spotcheck_min(@(x) -(x - 0.25).^3.*cosh(x), ...
    -0.75^3*cosh(1), pref);

%%
% Check operation for vectorized inputs.

fun_op = @(x) -[sin(10*x) airy(x) (x - 0.25).^3.*cosh(x)];
f = funcheb2(fun_op, pref);
[y, x] = min(f);
exact_max = -[1 airy(-1) 0.75^3*cosh(1)];
fx = -[sin(10*x(1)) airy(x(2)) (x(3) - 0.25).^3.*cosh(x(3))];
pass(6) = (all(abs(y - exact_max) < 10*f.epslevel) && ...
           all(abs(fx - exact_max) < 10*f.epslevel));

% Test for complex-valued funcheb2 objects.
pass(7) = test_spotcheck_min(@(x) exp(1i*x)-.5i*sin(x)+x, ...
    0.074968381369117 - 0.319744137826069i, pref);
    
fun_op = @(x) [exp(1i*x)-.5i*sin(x)+x (1+.5*(x-.1).^2).*exp(1i*x)];
f = funcheb2(fun_op, pref);
[y, x] = min(f);
exact_max = [0.074968381369117-0.319744137826069i, ...
    0.995004165278026+0.099833416646827i];
fx = fun_op(x); fx = fx([1 4]);
pass(8) = (all(abs(y - exact_max) < 100*f.epslevel) && ...
           all(abs(fx - exact_max) < 100*f.epslevel));

end

% Spot-check the results for a given function.
function result = test_spotcheck_min(fun_op, exact_min, pref)

f = funcheb2(fun_op, pref);
[y, x] = min(f);
fx = fun_op(x);
result = ((abs(y - exact_min) < 10*f.epslevel) && ... 
          (abs(fx - exact_min) < 10*f.epslevel));

end
