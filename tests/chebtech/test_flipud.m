% Test file for chebtech/flipud.

function pass = test_flipud(pref)

if ( nargin < 1 )
    pref = chebtech.pref;
end
tol = 10*pref.chebtech.eps;

pass = zeros(2, 4); % Pre-allocate pass matrix
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end
    
    % Try some standard calls to flipud 
    f = testclass.make(@(x) sin(x+.5), [], [], pref);
    g = testclass.make(@(x) sin(-x+.5), [], [], pref);
    h = flipud(f);
    pass(n, 1) = norm(g.values - h.values, inf) < tol;
    
    f = testclass.make(@(x) [sin(x+.5), exp(x)], [], [], pref);
    g = testclass.make(@(x) [sin(-x+.5), exp(-x)], [], [], pref);
    h = flipud(f);
    pass(n, 2) = norm(g.values - h.values, inf) < tol;
    
    f = testclass.make(@(x) sin(1i*x+.5), [], [], pref);
    g = testclass.make(@(x) sin(-1i*x+.5), [], [], pref);
    h = flipud(f);
    pass(n, 3) = norm(g.values - h.values, inf) < tol;
    
    f = testclass.make(@(x) [sin(x+.5), exp(1i*x)], [], [], pref);
    g = testclass.make(@(x) [sin(-x+.5), exp(-1i*x)], [], [], pref);
    h = flipud(f);
    pass(n, 4) = norm(g.values - h.values, inf) < tol;
end

end
