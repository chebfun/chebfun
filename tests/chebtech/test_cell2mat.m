% Test file for chebtech/cell2mat.m

function pass = test_cell2mat(pref)

if ( nargin < 2 )
    pref = chebtech.techPref();
end

for n = 1:2
    
    if (n == 1);
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], pref);
    g = testclass.make(@(x) sin(x), [], pref);
    h = testclass.make(@(x) [cos(x) exp(x)], [], pref);
    
    F = cell2mat([g h]);
    pass(n, 1) = all( sum(F - f) < max(f.vscale.*f.epslevel) );
end

end
