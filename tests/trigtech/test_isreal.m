% Test file for trigtech/isreal.m

function pass = test_isreal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Test a scalar-valued function.
f = testclass.make(@(x) sin(100*pi*x) + 1i*sin(cos(10*pi*x)), [], pref);
pass(1) = ~isreal(f);

f = testclass.make(@(x) 1i*cos(pi*x), [], pref);
pass(2) = ~isreal(f);

f = testclass.make(@(x) sin(pi*x), [], pref);
pass(3) = isreal(f);

% Test an array-valued function.
f = testclass.make(@(x) [sin(100*pi*x) + 1i*sin(cos(10*pi*x)), cos(10*pi*x)], [], pref);
pass(4) = ~isreal(f);

f = testclass.make(@(x) [1i*cos(pi*x), cos(sin(pi*x))], [], pref);
pass(5) = ~isreal(f);

f = testclass.make(@(x) [sin(20*pi*x), cos(sin(pi*x))], [], pref);
pass(6) = isreal(f);

end
