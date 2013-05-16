% Test file for fun/fliplr.m

function pass = test_fliplr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1,  2); % Pre-allocate pass matrix
for n = 1:1 % [TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();    
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    % Conduct a few very straightforward tests.
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    pass(n, 1) = isequal(f, fliplr(f));
    
    f = testclass.make(@(x) [sin(x) cos(x)], dom, [], [], pref);
    g = testclass.make(@(x) [cos(x) sin(x)], dom, [], [], pref);
    pass(n, 2) = isequal(fliplr(f), g);
end

end
