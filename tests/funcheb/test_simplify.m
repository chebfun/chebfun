function pass = test_simplify(pref)
% [TODO] Make this test more extensive.

% Get preferences:
if ( nargin < 1 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    % Set the tolerance:
    tol = 10*pref.funcheb.eps;
    
    %%
    % Test on a scalar-valued function:
    
    f = @(x) sin(x);
    pref.funcheb.n = 33;
    g = testclass.make(f, 0, pref);
    h = simplify(g);
    x = testclass.chebpts(14);
    pass(n, 1) = length(g) == 33 && length(h) == 14 && ...
        norm(f(x) - h.values, inf) < tol;
    
    %%
    % Test on a vector-valued function:
    
    f = @(x) [ sin(x), cos(x), exp(x) ];
    pref.funcheb.n = 33;
    g = testclass.make(f, 0, pref);
    h = simplify(g);
    x = testclass.chebpts(15);
    pass(n, 2) = length(g) == 33 && length(h) == 15 && ...
        norm(f(x) - h.values, inf) < tol;
end

end
