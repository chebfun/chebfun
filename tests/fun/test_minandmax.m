% Test file for fun/minandmax.m

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    %%
    % Spot-check the extrema for a few functions.
    pass(n, 1) = test_spotcheck_minmax(testclass, @(x) sin(10*x), dom, -1, 1, pref);
    pass(n, 2) = test_spotcheck_minmax(testclass, @(x) real(airy(x)), dom, 7.492128863997157e-07, 0.535656656015700, ...
        pref);
    pass(n, 3) = test_spotcheck_minmax(testclass, @(x) -1./(1 + x.^2), dom, -1, ...
        -0.02, pref);
    pass(n, 4) = test_spotcheck_minmax(testclass, ...
        @(x) (x/10).^3.*cosh(x/10), ...
        dom, (-.2)^3*cosh(-.2), 0.7^3*cosh(0.7), pref);

    %%
    % Check operation for array-valued inputs.
    fun_op = @(x) [sin(10*x) real(airy(x)) (x/10).^3.*cosh(x/10)];
    f = testclass.make(fun_op, dom, [], [], pref);
    [y, x] = minandmax(f);
    y_exact = [-1 7.492128863997157e-07  (-.2)^3*cosh(-.2);
                1 0.535656656015700 0.7^3*cosh(0.7)];
    pass(n, 5) = all(abs(y(:) - y_exact(:)) < 10*max(get(f, 'epslevel')));

    % Check that the points x are indeed extreme points of the function 
    % operator.
    pass(n, 6) = 1;
    for k = 1:1:size(f, 2)
        fx = fun_op(x(:, k));
        max(abs(fx(:, k)) - y_exact(:, k));
        if ( max(abs(fx(:, k) - y_exact(:, k))) > 10*get(f, 'epslevel') )
            pass(n, 6) = 0;
            break;
        end
    end

    %%  
    % Test complex-array-valued fun objects.
    f = testclass.make(@(x) [exp(sin(2*x)), 1i*cos(20*x)], dom);
    [vals, pos] = minandmax(f);
    f1 = testclass.make(@(x) exp(sin(2*x)), dom);
    [vals1, pos1] = minandmax(f1);
    f2 = testclass.make(@(x) 1i*cos(20*x), dom);
    [vals2, pos2] = minandmax(f2);
    pass(n, 7) = norm(abs(vals) - abs([vals1 vals2]), inf) < ...
        10*max(get(f, 'vscale').*get(f, 'epslevel'));
end

end

% Spot-check the results for a given function.
function result = test_spotcheck_minmax(testclass, fun_op, dom, exact_min, ...
    exact_max, pref)

f = testclass.make(fun_op, dom, [], [], pref);
[y, x] = minandmax(f);
y_exact = [exact_min ; exact_max];
fx = fun_op(x);
result = ((max(abs(y - y_exact)) < 10*get(f, 'epslevel')) && ... 
          (max(abs(fx - y_exact)) < 10*get(f, 'epslevel')));

end
