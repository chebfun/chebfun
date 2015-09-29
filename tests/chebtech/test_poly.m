% Test file for chebtech/poly.m

function pass = test_poly(pref)

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
    % Check a few simple examples.
    
    f = testclass.make(@(x) zeros(size(x)), [], pref);
    p = poly(f);
    pass(n, 1) = (norm(p, inf) <= 10*vscale(f)*eps);

    f = testclass.make(@(x) 3*ones(size(x)), [], pref);
    p = poly(f);
    pass(n, 2) = (norm(p - 3, inf) < 10*vscale(f)*eps);
    
    f = testclass.make(@(x) 6.4*x - 3i, [], pref);
    p = poly(f);
    pass(n, 3) = (norm(p - [6.4 (-3i)], inf) < 10*vscale(f)*eps);
    
    f = testclass.make(@(x) 2i*x.^5 - 3.2*x.^4 + 2*x.^2 - (1.2 + 3i), [], pref);
    p = poly(f);
    pass(n, 4) = (norm(p - [2i (-3.2) 0 2 0 -(1.2 + 3i)], inf) ...
        < 10*vscale(f)*eps);
    
    %%
    % Verify operation for array-valued chebtech objects.
    
    f = testclass.make(@(x) [3*ones(size(x)), (6.4*x - 3i), ... 
        (4*x.^2 - 2i*x + 3.7)], [], pref);
    p = poly(f);
    p_exact = [0 0     3;
               0 6.4   (-3i);
    	   4 (-2i) 3.7];
    pass(n, 5) = (norm(p(:) - p_exact(:), inf) < 10*max(vscale(f)*eps));
end

end
