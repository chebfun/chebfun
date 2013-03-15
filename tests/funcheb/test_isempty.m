function pass = test_isempty(varargin)

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    f = testclass.make();
    pass(n, 1) = isempty(f);
    
    f = testclass.make(@sin);
    pass(n, 2) = ~isempty(f);
    
    f = testclass.make(@(x) [sin(x), cos(x)]);
    pass(n, 3) = ~isempty(f);
    
    f = [ testclass.make(@sin), testclass.make(@sin) ];
    pass(n, 4) = ~isempty(f);
end

end
