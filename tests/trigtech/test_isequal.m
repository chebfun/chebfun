% Test file for trigtech/isequal.m

function pass = test_isequal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
f = testclass.make(@(x) sin(200*pi*x), [], pref);
g = f;
pass(1) = isequal(f, g) && isequal(g, f);

%%
g = testclass.make(@(x) cos(200*pi*x), [], pref);
pass(2) = ~isequal(f, g);

%%
g = testclass.make(@(x) [sin(200*pi*x), cos(200*pi*x)], [], pref);
pass(3) = ~isequal(f, g);

%%
f = g;
pass(4) = isequal(f, g);

%%
g = testclass.make(@(x) [sin(200*pi*x), exp(1i*200*pi*x)], [], pref);
pass(5) = ~isequal(f, g);

end
