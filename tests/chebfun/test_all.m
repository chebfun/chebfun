% Test file for @chebfun/all.m.

function pass = test_all(pref)

if (nargin < 1)
    pref = chebfunpref();
end

% Test scalar valued chebfuns.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
pass(1) = ~all(f);

f = chebfun(@(x) sin(x - 0.1), [-1 -0.5 0 0.5 1], pref);
pass(2) = ~all(f);

f = chebfun(@(x) exp(2*pi*1i*x), [-1 -0.5 0 0.5 1], pref);
pass(3) = all(f);

% Test array-valued chebfun.
f = chebfun(@(x) [sin(x) sin(x - 0.1) exp(2*pi*1i*x)], [-1 -0.5 0 0.5 1], pref);
pass(4) = isequal(all(f), logical([0 0 1]));

% Test on SINGFUN:
f = chebfun(@(x) sin(x)./(x+1), 'exps', [-1 0]);
pass(5) = ~all(f);

% Test for function defined on unbounded domain:

% Blowing-up functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

op = @(x) x.^2.*(1-exp(-x.^2))+3;
f = chebfun(op, dom, 'exps', [2 2]);
pass(6) = all(f);

end
