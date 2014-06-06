% Test file for chebtech/max.m

function pass = test_max(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    %%
    % Spot-check the extrema for a few functions.
    pass(n, 1) = test_spotcheck_max(testclass, @(x) ...
        ((x-0.2).^3 -(x-0.2) + 1).*sec(x-0.2), ...
        1.884217141925336, pref);
    pass(n, 2) = test_spotcheck_max(testclass, @(x) sin(10*x), 1, pref);
    pass(n, 3) = test_spotcheck_max(testclass, @airy, airy(-1), pref);
    pass(n, 4) = test_spotcheck_max(testclass, @(x) -1./(1 + x.^2), -0.5, pref);
    pass(n, 5) = test_spotcheck_max(testclass, @(x) (x - 0.25).^3.*cosh(x), ...
        0.75^3*cosh(1), pref);

    %%
    % Check operation for array-valued inputs.
    fun_op = @(x) [sin(10*x) airy(x) (x - 0.25).^3.*cosh(x)];
    f = testclass.make(fun_op, [], pref);
    [y, x] = max(f);
    exact_max = [1 airy(-1) 0.75^3*cosh(1)];
    fx = [sin(10*x(1)) airy(x(2)) (x(3) - 0.25).^3.*cosh(x(3))];
    pass(n, 6) = (all(abs(y - exact_max) < 10*f.epslevel) && ...
               all(abs(fx - exact_max) < 10*f.epslevel));

    %%
    % Test for complex-valued chebtech objects.
    pass(n, 7) = test_spotcheck_max(testclass, ...
        @(x) (x - 0.2).*(exp(1i*(x - 0.2))+1i*sin(x - 0.2)), ...
        -0.434829305372008 + 2.236893806321343i, pref);

    z = @(x) (x - 0.3 + 1i).^3 - 2i;
    fun_op = @(x) [sin(z(x)) sinh(z(x))];
    f = testclass.make(fun_op, [], pref);
    [y, x] = max(f);
    exact_max = [-10.017874927409903i 3.626860407847019];
    fx = [sin(z(x(1))) sinh(z(x(2)))];
    tol = 10*max(f.vscale, f.epslevel);
    pass(n, 8) = (all(abs(y - exact_max) < tol) && ...
                  all(abs(fx - exact_max) < tol));
end

end

% Spot-check the results for a given function.
function result = test_spotcheck_max(testclass, fun_op, exact_max, pref)

f = testclass.make(fun_op,[], pref);
[y, x] = max(f);
fx = fun_op(x);
result = (all(abs(y - exact_max) < 10*f.vscale.*f.epslevel) && ...
          all(abs(fx - exact_max) < 10*f.vscale.*f.epslevel));

end
