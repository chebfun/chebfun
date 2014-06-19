% Test file for bndfun/changeMap.m.

function pass = test_changeMap(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Test with two domains:
dom1 = [-3 1];
dom2 = [-2 7];
a = dom1(1);
b = dom1(2);
c = dom2(1);
d = dom2(2);

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom1) * rand(100, 1) + dom1(1);
y = c*(b - x)/(b - a) + d*(x - a)/(b - a);

f = bndfun(@(x) 1./(1 + x.^2), struct('domain', dom1), pref);

% Change from dom1 to dom2.
g = changeMap(f, dom2);
gv = feval(g, y);
fv = feval(f, x);
pass(1) = norm(gv - fv, inf) < 10*get(f, 'vscale')*get(f, 'epslevel');

% Change from dom2 to dom1.
f = changeMap(g, dom1);
fv = feval(f, x);
gv = feval(g, y);
pass(2) = norm(fv - gv, inf) < 10*get(f, 'vscale')*get(f, 'epslevel');

%% tests on singfun:

pref.blowup = 1;
op = @(x) 1./((x - a).*(x - b));
f = bndfun(op, struct('domain', dom1), pref);

% Change from dom1 to dom2.
g = changeMap(f, dom2);
gv = feval(g, y);
fv = feval(f, x);
pass(3) = all( abs(gv - fv) < 1e4*abs(op(x))*get(f, 'vscale')*get(f, 'epslevel') );

% Change from dom2 to dom1.
f = changeMap(g, dom1);
fv = feval(f, x);
gv = feval(g, y);
pass(4) = all( abs(fv - gv) < 1e4*abs(op(x))*get(f, 'vscale')*get(f, 'epslevel') );

end
