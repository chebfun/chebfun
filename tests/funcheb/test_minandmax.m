% Test file for funcheb/minandmax.

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    %%
    % Spot-check the extrema for a few functions.
    
    pass(n, 1) = test_spotcheck_minmax(testclass, @(x) ...
        ((x-0.2).^3 -(x-0.2) + 1).*sec(x-0.2), ...
        0.710869767377087, 1.884217141925336, pref);
    pass(n, 2) = test_spotcheck_minmax(testclass, @(x) sin(10*x), -1, 1, pref);
    pass(n, 3) = test_spotcheck_minmax(testclass, @airy, airy(1), airy(-1), ...
        pref);
    pass(n, 4) = test_spotcheck_minmax(testclass, @(x) -1./(1 + x.^2), -1, ...
        -0.5, pref);
    pass(n, 5) = test_spotcheck_minmax(testclass, ...
        @(x) (x - 0.25).^3.*cosh(x), ...
        (-1.25)^3*cosh(-1), 0.75^3*cosh(1), pref);

    %%
    % Check operation for vectorized inputs.
    
    fun_op = @(x) [sin(10*x) airy(x) (x - 0.25).^3.*cosh(x)];
    f = testclass.make(fun_op, pref);
    [y, x] = minandmax(f);
    y_exact = [-1 airy(1)  (-1.25)^3*cosh(-1);
                1 airy(-1) 0.75^3*cosh(1)];

    pass(n, 6) = all(abs(y(:) - y_exact(:)) < 10*f.epslevel);

    % Check that the points x are indeed extreme points of the function 
    % operator.
    for k = 1:1:size(f.coeffs, 2)
        fx = fun_op(x(:, k));
        if ( max(abs(fx(:, k) - y_exact(:, k))) > 10*f.epslevel )
            pass(n, 6) = 0;
            break;
        end
    end
end

end

% Spot-check the results for a given function.
function result = test_spotcheck_minmax(testclass, fun_op, exact_min, ...
    exact_max, pref)

f = testclass.make(fun_op, pref);
[y, x] = minandmax(f);
y_exact = [exact_min ; exact_max];
fx = fun_op(x);
result = ((max(abs(y - y_exact)) < 10*f.epslevel) && ... 
          (max(abs(fx - y_exact)) < 10*f.epslevel));

end
