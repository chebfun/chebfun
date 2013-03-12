function pass = test_mat2cell(pref)

if ( nargin < 2 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], pref);
    g = testclass.make(@(x) sin(x), pref);
    h = testclass.make(@(x) [cos(x) exp(x)], pref);
    
    F = mat2cell(f, 1, [1 2]);
    pass(n, 1) = sum(F(1) - g) < g.epslevel;
    pass(n, 2) = all( sum(F(1) - g) < h.vscale*h.epslevel );
    % [TODO]:  This test fails for funcheb1.
    pass(n, 3) = abs(F(1).vscale - sin(1)) < 10*F(1).epslevel;
    % [TODO]:  This test fails for funcheb1.
    pass(n, 4) = all( abs(F(2).vscale - [cos(0) exp(1)]) < 10*F(2).epslevel );
end

end
