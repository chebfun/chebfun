% Test file for funcheb2/minandmax.

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

%%
% Spot-check the extrema for a few functions.

pass(1) = test_spotcheck_minmax(@(x) ((x-0.2).^3 -(x-0.2) + 1).*sec(x-0.2), ...
    0.710869767377087, 1.884217141925336, pref);
pass(2) = test_spotcheck_minmax(@(x) sin(10*x), -1, 1, pref);
pass(3) = test_spotcheck_minmax(@airy, airy(1), airy(-1), pref);
pass(4) = test_spotcheck_minmax(@(x) -1./(1 + x.^2), -1, -0.5, pref);
pass(5) = test_spotcheck_minmax(@(x) (x - 0.25).^3.*cosh(x), ...
    (-1.25)^3*cosh(-1), 0.75^3*cosh(1), pref);

%%
% Check operation for vectorized inputs.

fun_op = @(x) [sin(10*x) airy(x) (x - 0.25).^3.*cosh(x)];
f = funcheb2(fun_op, pref);
[y, x] = minandmax(f);
y_exact = [-1 airy(1)  (-1.25)^3*cosh(-1);
            1 airy(-1) 0.75^3*cosh(1)];

pass(6) = all(abs(y(:) - y_exact(:)) < 10*f.epslevel);

% Check that the points x are indeed extreme points of the function operator.
for k = 1:1:size(f.coeffs, 2)
    fx = fun_op(x(:, k));
    if ( max(abs(fx(:, k) - y_exact(:, k))) > 10*f.epslevel )
        pass(6) = 0;
        break;
    end
end

% Test for complex-valued funcheb2 objects:
% These are tested sufficiently in test_min and test_max.


end

% Spot-check the results for a given function.
function result = test_spotcheck_minmax(fun_op, exact_min, exact_max, pref)

f = funcheb2(fun_op, pref);
[y, x] = minandmax(f);
y_exact = [exact_min ; exact_max];
fx = fun_op(x);
result = ((max(abs(y - y_exact)) < 10*f.epslevel) && ... 
          (max(abs(fx - y_exact)) < 10*f.epslevel));

end
