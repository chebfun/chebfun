% Test file for chebtech/simplify.m

function pass = test_simplify(pref)

% SIMPLIFY() is tested quite extensively in TEST_HAPPINESSCHECK() and
% TEST_PROLONG().

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.pref;
end

pass = zeros(2, 6); % Pre-allocate pass matrix
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Set the tolerance:
    tol = 10*pref.chebtech.eps;
    
    %%
    % Test on a scalar-valued function:
    
    f = @(x) sin(x);
    pref.chebtech.n = 33;
    g = testclass.make(f, [], [], pref);
    h = simplify(g);
    x = testclass.chebpts(length(h));
    pass(n, 1) = length(g) == 33;
    pass(n, 2) = length(h) == 14;
    pass(n, 3) = norm(f(x) - h.values, inf) < tol;
    
    %%
    % Test on a array-valued function:
    
    f = @(x) [ sin(x), cos(x), exp(x) ];
    pref.chebtech.n = 33;
    g = testclass.make(f, [], [], pref);
    h = simplify(g);
    x = testclass.chebpts(length(h));
    pass(n, 4) = length(g) == 33;
    pass(n, 5) = abs(length(h) - 15) < 2;
    pass(n, 6) = norm(f(x) - h.values, inf) < tol;
end

end
