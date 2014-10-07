% Test file for trigtech/isnan.m

function pass = test_isnan(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Test a scalar-valued function.
f = testclass.make(@(x) cos(pi*x), [], pref);
pass(1) = ~isnan(f);

% Test an array-valued function.
f = testclass.make(@(x) [cos(pi*x), cos(pi*x).^2], [], pref);
pass(2) = ~isnan(f);

% Artificially construct and test a NaN-valued function.
f = testclass.make(NaN);
pass(3) = isnan(f);

% Test a NaN scalar-valued function.
try
    f = testclass.make(@(x) cos(pi*x) + NaN);
    pass(4) = isnan(f);
catch ME
    pass(4) = strcmpi(ME.message, 'Cannot handle functions that evaluate to Inf or NaN.');
end

% Test a NaN array-valued function.
try
    f = testclass.make(@(x) [cos(pi*x) + NaN, cos(pi*x)]);
    pass(5) = isnan(f);
catch ME
    pass(5) = strcmpi(ME.message, 'Cannot handle functions that evaluate to Inf or NaN.');
end

end
