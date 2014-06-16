% Test file for chebtech/flipud.m

function pass = test_flipud(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end
    
    % Try some standard calls to flipud 
    f = testclass.make(@(x) sin(x+.5), [], pref);
    g = testclass.make(@(x) sin(-x+.5), [], pref);
    h = flipud(f);
    pass(n, 1) = norm(g.coeffs - h.coeffs, inf) < 10*h.vscale.*h.epslevel;
    
    f = testclass.make(@(x) [sin(x+.5), exp(x)], [], pref);
    g = testclass.make(@(x) [sin(-x+.5), exp(-x)], [], pref);
    h = flipud(f);
    pass(n, 2) = norm(g.coeffs - h.coeffs, inf) < 10*max(h.vscale.*h.epslevel);
    
    f = testclass.make(@(x) sin(1i*x+.5), [], pref);
    g = testclass.make(@(x) sin(-1i*x+.5), [], pref);
    h = flipud(f);
    pass(n, 3) = norm(g.coeffs - h.coeffs, inf) < 10*h.vscale.*h.epslevel;
    
    f = testclass.make(@(x) [sin(x+.5), exp(1i*x)], [], pref);
    g = testclass.make(@(x) [sin(-x+.5), exp(-1i*x)], [], pref);
    h = flipud(f);
    pass(n, 4) = norm(g.coeffs - h.coeffs, inf) < 10*max(h.vscale.*h.epslevel);
end

end
