% Test file for @chebtech/assignColumns.m.

function pass = test_assignColumns(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], pref);

    g = testclass.make(@(x) [x x.^2], [], pref);
    h = assignColumns(f, [1 3], g);
    h_exact = @(x) [x cos(x) x.^2];
    err = feval(h, x) - h_exact(x);
    pass(1) = h.ishappy && (norm(err(:), inf) < 10*max(vscale(h)*eps));

    %
    g = testclass.make(@(x) sqrt(x), [], pref);
    h = assignColumns(f, 1, g);
    pass(2) = ~h.ishappy;
    
    %
    h = assignColumns(f, 1, []);
    pass(3) = all(size(vscale(h)) == [1, 2]);
    
end

end
