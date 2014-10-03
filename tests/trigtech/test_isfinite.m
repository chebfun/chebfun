% Test file for fourtech/isfinite.m

function pass = test_isfinite(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fourtech.techPref();
end

testclass = fourtech();

%%
% Test a scalar-valued function.
y = ones(11, 1);
y(4) = inf;
f = testclass.make({[], y}, [], []);
pass(1) = ~isfinite(f);

%%
% Test an array-valued function.
y = ones(11,1);
y(4) = inf;
f = testclass.make({[], y}, [], pref);
pass(2) = ~isfinite(f);

%%
% Test a finite scalar-valued function.
p = pref;
f = testclass.make(@(x) x, [], pref);
pass(3) = isfinite(f);

%%
% Test a finite array-valued function.
f = testclass.make(@(x) [x, x], [], pref);
pass(4) = isfinite(f);

end
