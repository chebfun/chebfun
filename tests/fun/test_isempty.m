% Test file for fun/isempty.m

function pass = test_isempty(varargin)

pass = zeros(1, 4); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    f = testclass.make();
    pass(n, 1) = isempty(f);
    
    f = testclass.make(@sin, dom);
    pass(n, 2) = ~isempty(f);
    
    f = testclass.make(@(x) [sin(x), cos(x)], dom);
    pass(n, 3) = ~isempty(f);
    
    f = [ testclass.make(@sin, dom), testclass.make(@sin, dom) ];
    pass(n, 4) = ~isempty(f);
end

end