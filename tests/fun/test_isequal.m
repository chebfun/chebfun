% Test file for fun/isequal.m

function pass = test_isequal(pref)

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
        dom = [0 inf];        
    end

    % Run a few very straightforward tests.
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    g = f;
    pass(n, 1) = isequal(f, g) && isequal(g, f);
    
    g = testclass.make(@(x) cos(x), dom, [], [], pref);
    pass(n, 2) = ~isequal(f, g);
    
    g = testclass.make(@(x) [sin(x), cos(x)], dom, [], [], pref);
    pass(n, 3) = ~isequal(f, g);
    
    f = g;
    pass(n, 4) = isequal(f, g);
    
    g = testclass.make(@(x) [sin(x), exp(x)], dom, [], [], pref);
    pass(n, 5) = ~isequal(f, g);
    
    %% Test on singular function:
    
    pow1 = -0.5;
    pow2 = -0.6;
    op1 = @(x) (x - dom(2)).^pow1.*sin(x);
    op2 = @(x) (x - dom(2)).^pow2.*(cos(x).^2+1);
    pref.singPrefs.exponents = [0 pow1];
    f = bndfun(op1, dom, [], [], pref);
    pref.singPrefs.exponents = [0 pow2];
    g = bndfun(op2, dom, [], [], pref);
    pass(n, 6) = ~isequal(f, g);
    
end

end
