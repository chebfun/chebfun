% Test file for funcheb/mat2cell.m

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

    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], [], pref);
    g = testclass.make(@(x) sin(x), [], [], pref);
    h = testclass.make(@(x) [cos(x) exp(x)], [], [], pref);
    
    F = mat2cell(f, 1, [1 2]);
    pass(n, 1) = sum(F(1) - g) < g.epslevel;
    pass(n, 2) = all( sum(F(1) - g) < h.vscale*h.epslevel );
end

end
