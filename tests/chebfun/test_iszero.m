function pass = test_iszero(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Test scalars:
f = chebfun(0, pref);
pass(1) = iszero(f);

f = chebfun([], pref);
pass(2) = iszero(f);

f = chebfun(2, pref);
pass(3) = ~iszero(f);

% Test piecewise domains:
f = chebfun(0, [-1, 0, 1], pref);
pass(4) = iszero(f);

f = chebfun([], [-1, 0, 1], pref);
pass(5) = iszero(f);

f = chebfun({0, 1}, [-1, 0, 1], pref);
pass(6) = ~iszero(f);

f = chebfun({1, 0}, [-1, 0, 1], pref);
pass(7) = ~iszero(f);

% Test arrays:
f = chebfun([0 0], pref);
pass(8) = all(iszero(f));

f = chebfun([0 1], pref);
pass(9) = all(iszero(f) == [1 0]);

%% Test on singular function:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) sin(100*x).*(x-dom(1)).^pow, dom, 'exps', [pow 0], ...
    'splitting', 'on');
pass(10) = ~iszero(f);

%% Test for functions defined on unbounded domain:

% Set the domain:
dom = [-Inf Inf];

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
f = chebfun(op, dom, 'exps', [2 2]);
pass(11) = ~iszero(f);

% Function defined on [0 Inf]:

% Specify the domain:
dom = [0 Inf];

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
pass(12) = ~iszero(f);

% Function defined on [0 Inf]:

% Set the domain:
dom = [-Inf -3*pi];

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = chebfun(op, dom, 'exps', [0 -1]); 
pass(13) = ~iszero(f);

% Zero function:
f = chebfun(@(x) 0*x, dom); 
pass(14) = iszero(f);

% [TODO]: Add these once SUBSREF is implemented.
% f = chebfun(0, pref);
% f(0) = 1;
% pass(10) = iszero(f);
% 
% f = chebfun([0 0], pref);
% f(0) = [0, 1];
% pass(11) = all(iszero(f) == [1, 0]);

end
