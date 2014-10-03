% Test file for trigtech/any.m.

function pass = test_any(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Check behavior for any() down columns.
pass(1) = ~any(testclass());

f = testclass.make(@(x) [sin(pi*x) 0*x cos(pi*x)]);
pass(2) = isequal(any(f), [1 0 1]);

% Check behavior for any() across rows.
g = any(f, 2);
pass(3) = isequal(g.coeffs, 1);

f = testclass.make(@(x) [0*x 0*x]);
g = any(f, 2);
pass(4) = isequal(g.coeffs, 0);

end
