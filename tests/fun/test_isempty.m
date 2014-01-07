% Test file for fun/isempty.m

function pass = test_isempty(varargin)

for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
        dom = [0 inf];
    end

    % Test an empty FUN:
    f = testclass.make();
    pass(n, 1) = isempty(f);
    
    % Test an non-empty FUN:
    f = testclass.make(@sin, dom);
    pass(n, 2) = ~isempty(f);
    
    % Test an non-empty array-valued FUN:
    f = testclass.make(@(x) [sin(x), cos(x)], dom);
    pass(n, 3) = ~isempty(f);
    
    % Test an array of FUN objects:
    f = [ testclass.make(@sin, dom), testclass.make(@sin, dom) ];
    pass(n, 4) = ~isempty(f);

    % Test an array of empty FUN objects:
    f = [ testclass.make(), testclass.make() ];
    pass(n, 5) = isempty(f);
    
    %% Integration of singfun:
    
    pow = -0.5;
    op = @(x) (x - dom(2)).^pow.*sin(x);
    pref.singPrefs.exponents = [0 pow];
    f = testclass.make(op, dom, [], [], pref);
    pass(n, 6) = ~isempty(f);
    
end

end
