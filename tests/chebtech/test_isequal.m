% Test file for chebtech/isequal.m

function pass = test_isequal(pref)

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
    % Run a few very straightforward tests.
    
    f = testclass.make(@(x) sin(x), [], pref);
    g = f;
    pass(n, 1) = isequal(f, g) && isequal(g, f);
    
    g = testclass.make(@(x) cos(x), [], pref);
    pass(n, 2) = ~isequal(f, g);
    
    g = testclass.make(@(x) [sin(x) cos(x)], [], pref);
    pass(n, 3) = ~isequal(f, g);
    
    f = g;
    pass(n, 4) = isequal(f, g);
    
    g = testclass.make(@(x) [sin(x) exp(x)], [], pref);
    pass(n, 5) = ~isequal(f, g);
end

end
