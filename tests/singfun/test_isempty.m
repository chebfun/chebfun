% Test file for singfun/isempty.m

function pass = test_isempty(pref)

% No preferences used.

% create an empty SINGFUN
f = singfun;

% create a zero SINGFUN
g = singfun.zeroSingFun();

% create a non-zero SINGFUN
h = singfun(@(x)1./(1+x), [-1, 0] );

% Test
pass(1) = isempty(f);
pass(2) = ~isempty(g);
pass(3) = ~isempty(h);
