% Test file for trigtech/isinf.m

function pass = test_isinf(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Test a scalar-valued function.
y = ones(11,1);
y(4) = inf;
f = testclass.make({[],y}, [], []);
pass(1) = isinf(f);

% Test an array-valued function.
y = ones(11,1);
y(4) = inf;
f = testclass.make({[],y}, [], []);
pass(2) = isinf(f);

% Test a finite scalar-valued function.
f = testclass.make(@(x) cos(pi*x), [], pref);
pass(3) = ~isinf(f);

% Test a finite array-valued function.
f = testclass.make(@(x) [cos(pi*x), cos(pi*x)], [], pref);
pass(4) = ~isinf(f);

end
