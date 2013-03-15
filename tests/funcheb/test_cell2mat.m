function pass = test_cell2mat(pref)

if ( nargin < 2 )
    pref = funcheb.pref;
end

for (n = 1:2)
    if (n == 1);
        testclass = funcheb1();
    else
        testclass = funcheb2();
    end

    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], [], pref);
    g = testclass.make(@(x) sin(x), [], [], pref);
    h = testclass.make(@(x) [cos(x) exp(x)], [], [], pref);
    
    F = cell2mat([g h]);
    pass(n, 1) = all( sum(F - f) < f.vscale*f.epslevel );
end

end
