% Test file for chebtech/fliplr.m

function pass = test_fliplr(pref)

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
    % Conduct a few very straightforward tests.
    
    f = testclass.make(@(x) sin(x), [], pref);
    pass(n, 1) = isequal(f, fliplr(f));
    
    f = testclass.make(@(x) [sin(x) cos(x)], [], pref);
    g = testclass.make(@(x) [cos(x) sin(x)], [], pref);
    pass(n, 2) = isequal(fliplr(f), g);
end

end
