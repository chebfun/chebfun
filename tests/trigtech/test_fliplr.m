% Test file for trigtech/fliplr.m

function pass = test_fliplr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
% Conduct a few very straightforward tests.

f = testclass.make(@(x) cos(pi*x), [], pref);
pass(1) = isequal(f, fliplr(f));

f = testclass.make(@(x) [sin(pi*x) cos(pi*x)], [], pref);
g = testclass.make(@(x) [cos(pi*x) sin(pi*x)], [], pref);
pass(2) = isequal(fliplr(f), g);

end
