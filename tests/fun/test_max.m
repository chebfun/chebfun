% Test file for fun/max.m

function pass = test_max(pref)

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
    pass(n, 1) = test_spotcheck_max(testclass, @(x) sin(10*x), dom, 1, pref);
    pass(n, 2) = test_spotcheck_max(testclass, @airy, dom, 0.535656656015700, pref);
    pass(n, 3) = test_spotcheck_max(testclass, @(x) -1./(1 + x.^2), dom, -.02, pref);
    pass(n, 4) = test_spotcheck_max(testclass, @(x) (x/10).^3.*cosh(x/10), ...
        dom, 0.7^3*cosh(0.7), pref);

    %%
    % Check operation for array-valued inputs.
    fun_op = @(x) [sin(10*x) airy(x) (x/10).^3.*cosh(x/10)];
    f = testclass.make(fun_op, dom, [], [], pref);
    [y, x] = max(f);
    exact_max = [1 0.535656656015700 0.7^3*cosh(0.7)];
    fx = [sin(10*x(1)) airy(x(2)) (x(3)/10).^3.*cosh(x(3)/10)];
    pass(n, 5) = (all(abs(y - exact_max) < 10*f.onefun.epslevel) && ...
               all(abs(fx - exact_max) < 10*f.onefun.epslevel));

    %%
    % Test for complex-valued fun objects.
    pass(n, 6) = test_spotcheck_max(testclass, ...
        @(x) (x/2).*(exp(1i*(x/2))+1i*sin(x/2)), dom, ...
         -3.277598405517787 - 2.455482593827339i, pref);

    fun_op = @(x) [((x-2).^2/4+1).*exp(1i*(x/2)) ... 
                  -((x+1).^2/4+1).*exp(1i*(x/2))];
    f = testclass.make(fun_op, dom, [], [], pref);
    [y, x] = max(f);
    exact_max = [-6.789310982858273-2.543178400749744i 15.919763683943538+5.963314870723537i];
    fx = [((x(1)-2).^2/4+1).*exp(1i*(x(1)/2)) ... 
         -((x(2)+1).^2/4+1).*exp(1i*(x(2)/2))];
    pass(n, 7) = (all(abs(y - exact_max) < 100*f.onefun.epslevel) && ...
                  all(abs(fx - exact_max) < 100*f.onefun.epslevel));
end

end

% Spot-check the results for a given function.
function result = test_spotcheck_max(testclass, fun_op, dom, exact_max, pref)

f = testclass.make(fun_op, dom, [], [], pref);
[y, x] = max(f);
fx = fun_op(x);
result = (all(abs(y - exact_max) < 10*f.onefun.epslevel) && ...
          all(abs(fx - exact_max) < 10*f.onefun.epslevel));

end
