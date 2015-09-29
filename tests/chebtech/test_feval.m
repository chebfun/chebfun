% Test file for chebtech/feval.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(7681);
x = 2 * rand(1000, 1) - 1;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Spot-check values for a couple of functions.  We can only expect 
    % accuracy on % the order of the truncation level, so we use this as our 
    % criterion.
    
    f = testclass.make(@(x) exp(x) - 1, [], pref);
    f_exact = @(x) exp(x) - 1;
    pass(n, 1) = (norm(feval(f, x) - f_exact(x), inf) < ...
        10*vscale(f)*eps);
    
    f = testclass.make(@(x) 1./(1 + x.^2), [], pref);
    f_exact = @(x) 1./(1 + x.^2);
    pass(n, 2) = (norm(feval(f, x) - f_exact(x), inf) < ...
        10*vscale(f)*eps);
    
    f = testclass.make(@(x) cos(1e4*x), [], pref);
    f_exact = @(x) cos(1e4*x);
    pass(n, 3) = (norm(feval(f, x) - f_exact(x), inf) < ...
        1e4*vscale(f)*eps);
        
    
    z = exp(2*pi*1i/6);
    f = testclass.make(@(t) sinh(t*z), [], pref);
    f_exact = @(t) sinh(t*z);
    pass(n, 4) = (norm(feval(f, x) - f_exact(x), inf) < ...
        10*vscale(f)*eps);
    
    %%
    % Check row vector and matrix input.
    
    err = feval(f, x.') - f_exact(x.');
    pass(n, 5) = (all(size(err) == [1 1000])) && (norm(err(:), inf) < ...
        10*vscale(f)*eps);
    
    x_mtx = reshape(x, [100 10]);
    err = feval(f, x_mtx) - f_exact(x_mtx);
    pass(n, 6) = (all(size(err) == [100 10])) && (norm(err(:), inf) < ...
        10*vscale(f)*eps);
    
    x_3mtx = reshape(x, [10 10 10]);
    err = feval(f, x_3mtx) - f_exact(x_3mtx);
    pass(n, 7) = (all(size(err) == [10 10 10])) && (norm(err(:), inf) < ...
        10*vscale(f)*eps);
    
    %%
    % Check operation for array-valued chebtech objects.
    
    f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], pref);
    f_exact = @(x) [sin(x) x.^2 exp(1i*x)];
    err = feval(f, x) - f_exact(x);
    pass(n, 8) = all(max(abs(err)) < 10*max(vscale(f)*eps));
    
    %%
    % Test for evaluating array-valued chebtech objects at matrix arguments if 
    % the operation makes sense.
    f = testclass.make(@(x) [sin(pi*x) cos(pi*x) exp(pi*x)], [], pref);
    x2 = [-1 0 1 ; .25 .5 .75];
    fx = feval(f, x2);
    f_exact = [0 0 0 -1 1 -1 exp(-pi) 1 exp(pi)
              [1 sqrt(2) 1 1 0 -1]/sqrt(2) exp(pi.*[.25 .5 .75])];
    pass(n, 9) = all(all(abs(fx - f_exact) < 10*max(vscale(f)*eps)));
end

end
