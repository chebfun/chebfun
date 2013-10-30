% Test file for bndfun/changeMap.m.

function pass = test_changeMap(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebpref();
end

% Test with two domains:
dom1 = [-3 1];
dom2 = [-2 7];
a = dom1(1);
b = dom1(2);
c = dom2(1);
d = dom2(2);

f = bndfun(@(x) 1./(1 + x.^2), dom1, [], [], pref);

% Change from dom1 to dom2.
g = changeMap(f, dom2);
y = ((d - c)/2) * rand(100, 1) + (d + c)/2;
gv = feval(g, y);
x = a*(d - y)/(d - c) + b*(y - c)/(d - c);
fv = feval(f, x);
pass(1) = norm(gv - fv, inf) < 10*get(f, 'vscale')*get(f, 'epslevel');

% Change from dom2 to dom1.
f = changeMap(g, dom1);
y = ((b - a)/2) * rand(100, 1) + (b + a)/2;
fv = feval(f, y);
x = c*(b - y)/(b - a) + d*(y - a)/(b - a);
gv = feval(g, x);
pass(2) = norm(fv - gv, inf) < 10*get(f, 'vscale')*get(f, 'epslevel');

end
