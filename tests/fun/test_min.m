% Test file for fun/min.m

function pass = test_min(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 7); % Pre-allocate pass matrix.
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else
        testclass = unbndfun();
    end
    
    %%
    % Spot-check the extrema for a few functions.
    pass(n, 1) = test_spotcheck_min(testclass, @(x) -sin(10*x), dom, -1, pref);
    pass(n, 2) = test_spotcheck_min(testclass, @(x) -airy(x), dom, -0.535656656015700, pref);
    pass(n, 3) = test_spotcheck_min(testclass, @(x) 1./(1 + x.^2), dom, 0.02, pref);
    pass(n, 4) = test_spotcheck_min(testclass, @(x) -(x/10).^3.*cosh(x/10), ...
        dom, -0.7^3*cosh(0.7), pref);
    
    %%
    % Check operation for vectorized inputs.
    
    fun_op = @(x) -[sin(10*x) airy(x) (x/10).^3.*cosh(x/10)];
    f = testclass.make(fun_op, dom, [], [], pref);
    [y, x] = min(f);
    exact_max = -[1 0.535656656015700 0.7^3*cosh(0.7)];
    fx = -[sin(10*x(1)) airy(x(2)) (x(3)/10).^3.*cosh(x(3)/10)];
    pass(n, 5) = (all(abs(y - exact_max) < 10*f.onefun.epslevel) && ...
               all(abs(fx - exact_max) < 10*f.onefun.epslevel));
    
    %%
    % Test for complex-valued chebtech objects.
    pass(n, 6) = test_spotcheck_min(testclass, ...
        @(x) (x/2).*(exp(1i*(x/2))+1i*sin(x/2)), dom, ...
        0, pref);
        
    fun_op = @(x) [((x-2).^2/4+1).*exp(1i*(x/2)) ... 
                  -((x+1).^2/4+1).*exp(1i*(x/2))];
    f = testclass.make(fun_op, dom, [], [], pref);
    [y, x] = min(f);
    exact_max = [exp(1i) -exp(-1i/2)];
    fx = fun_op(x); 
    fx = fx([1 4]);
    pass(n, 7) = (all(abs(y - exact_max) < 100*f.onefun.epslevel) && ...
                  all(abs(fx - exact_max) < 100*f.onefun.epslevel));
end

end

% Spot-check the results for a given function.
function result = test_spotcheck_min(testclass, fun_op, dom, exact_min, pref)

f = testclass.make(fun_op, dom, [], [], pref);
[y, x] = min(f);
fx = fun_op(x);
result = ((abs(y - exact_min) < 100*f.onefun.epslevel) && ... 
          (abs(fx - exact_min) < 100*f.onefun.epslevel));

end
