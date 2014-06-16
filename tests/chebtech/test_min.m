% Test file for chebtech/min.m

function pass = test_min(pref)

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
    
    pass(n, 1) = test_spotcheck_min(testclass, @(x) ...
        -((x-0.2).^3 -(x-0.2) + 1).*sec(x-0.2), ...
        -1.884217141925336, pref);
    pass(n, 2) = test_spotcheck_min(testclass, @(x) -sin(10*x), -1, pref);
    pass(n, 3) = test_spotcheck_min(testclass, @(x) -airy(x), -airy(-1), pref);
    pass(n, 4) = test_spotcheck_min(testclass, @(x) 1./(1 + x.^2), 0.5, pref);
    pass(n, 5) = test_spotcheck_min(testclass, @(x) -(x - 0.25).^3.*cosh(x), ...
        -0.75^3*cosh(1), pref);
    
    %%
    % Check operation for array-valued inputs.
    
    fun_op = @(x) -[sin(10*x) airy(x) (x - 0.25).^3.*cosh(x)];
    f = testclass.make(fun_op, [], pref);
    [y, x] = min(f);
    exact_max = -[1 airy(-1) 0.75^3*cosh(1)];
    fx = -[sin(10*x(1)) airy(x(2)) (x(3) - 0.25).^3.*cosh(x(3))];
    tol = 10*max(f.vscale.*f.epslevel);
    pass(n, 6) = (all(abs(y - exact_max) < tol) && ...
               all(abs(fx - exact_max) < tol));
    
    % Test for complex-valued chebtech objects.
    pass(n, 7) = test_spotcheck_min(testclass, @(x) exp(1i*x)-.5i*sin(x)+x, ...
        0.074968381369117 - 0.319744137826069i, pref);
        
    fun_op = @(x) [exp(1i*x)-.5i*sin(x)+x (1+.5*(x-.1).^2).*exp(1i*x)];
    f = testclass.make(fun_op, [], pref);
    [y, x] = min(f);
    exact_max = [0.074968381369117-0.319744137826069i, ...
        0.995004165278026+0.099833416646827i];
    fx = fun_op(x); fx = fx([1 4]);
    tol = 10*max(f.vscale.*f.epslevel);
    pass(n, 8) = (all(abs(y - exact_max) < tol) && ...
                  all(abs(fx - exact_max) < tol));
end

end

% Spot-check the results for a given function.
function result = test_spotcheck_min(testclass, fun_op, exact_min, pref)

f = testclass.make(fun_op, [], pref);
[y, x] = min(f);
fx = fun_op(x);
result = ((abs(y - exact_min) < 10*f.vscale.*f.epslevel) && ... 
          (abs(fx - exact_min) < 10*f.vscale.*f.epslevel));

end
