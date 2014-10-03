% Test file for @trigtech/assignColumns.m.

function pass = test_assignColumns(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

testclass = trigtech();

f = testclass.make(@(x) [sin(pi*x) cos(pi*x) exp(1i*pi*x)], [], pref);

g = testclass.make(@(x) [exp(cos(pi*x)) exp(sin(pi*x))], [], pref);
h = assignColumns(f, [1 3], g);
h_exact = @(x) [exp(cos(pi*x)) cos(pi*x) exp(sin(pi*x))];
err = feval(h, x) - h_exact(x);
pass(1) = h.ishappy && (norm(err(:), inf) < 3e2*max(h.vscale.*h.epslevel));

g = testclass.make(@(x) x, [], pref);
h = assignColumns(f, 1, g);
pass(2) = ~h.ishappy;

end
