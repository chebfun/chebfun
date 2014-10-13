% Test file for trigtech/size.m

function pass = test_size(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

f = testclass.make(@(x) sin(10*pi*x), [], pref);
pass(1) = all(size(f) == size(f.coeffs));

f = testclass.make(@(x) [sin(19*pi*x), cos(22*pi*x), exp(10i*pi*x)], [], pref);
pass(2) = all(size(f) == size(f.coeffs));

p = pref;
p.fixedLength = 101;
f = testclass.make(@(x) [sin(19*pi*x), cos(22*pi*x), exp(10i*pi*x)], [], p);
pass(3) = all(size(f) == [101, 3]);

end
