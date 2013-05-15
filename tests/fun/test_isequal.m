% Test file for fun/isequal.m

function pass = test_isequal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 5); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
    else 
        testclass = unbndfun();
    end

    %%
    % Run a few very straightforward tests.
    
    dom = [-2 7];
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    g = f;
    pass(n, 1) = isequal(f, g) && isequal(g, f);
    
    g = testclass.make(@(x) cos(x), dom, [], [], pref);
    pass(n, 2) = ~isequal(f, g);
    
    g = testclass.make(@(x) [sin(x) cos(x)], dom, [], [], pref);
    pass(n, 3) = ~isequal(f, g);
    
    f = g;
    pass(n, 4) = isequal(f, g);
    
    g = testclass.make(@(x) [sin(x) exp(x)], dom, [], [], pref);
    pass(n, 5) = ~isequal(f, g);
end

end
