% Test file for funcheb/fliplr.

function pass = test_fliplr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    %%
    % Conduct a few very straightforward tests.
    
    f = testclass.make(@(x) sin(x), [], [], pref);
    pass(n, 1) = isequal(f, fliplr(f));
    
    f = testclass.make(@(x) [sin(x) cos(x)], [], [], pref);
    g = testclass.make(@(x) [cos(x) sin(x)], [], [], pref);
    pass(n, 2) = isequal(fliplr(f), g);
end

end
