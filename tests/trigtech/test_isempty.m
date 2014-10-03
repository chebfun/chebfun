% Test file for trigtech/isempty.m

function pass = test_isempty(varargin)

testclass = trigtech();

%%
f = testclass.make();
pass(1) = isempty(f);

%%
f = testclass.make(@(x) sin(200*pi*x));
pass(2) = ~isempty(f);

%%
f = testclass.make(@(x) [sin(200*pi*x), cos(200*pi*x)]);
pass(3) = ~isempty(f);

%%
f = [ testclass.make(@(x) sin(200*pi*x)), testclass.make(@(x) sin(200*pi*x)) ];
pass(4) = ~isempty(f);

%%
f = [ testclass.make() testclass.make() ];
pass(5) = isempty(f);

end