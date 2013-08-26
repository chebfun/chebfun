% Test file for fun/isempty.m

function pass = test_isempty(varargin)

pass = zeros(1, 4); % Pre-allocate pass matrix
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
end

end