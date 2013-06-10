% Test file for bndfun constructor.

function pass = test_changeMap(pref)

% Get preferences:
if ( nargin < 1 )
    pref = bndfun.pref;
end
% Set the tolerance:
tol = pref.bndfun.eps;

% test with two domains
pass = zeros(1, 2); % Pre-allocate pass matrix.
dom1 = [-3 1];
dom2 = [-2 7];
a = dom1(1);
b = dom1(2);
c = dom2(1);
d = dom2(2);

% forth
f = bndfun(@(x) 1./(1+x.^2), dom1, [], [], pref);
g = changeMap(f, dom2);
y = ((d-c)/2) * rand(100, 1) + (d+c)/2;
gv = feval(g,y);
x = a*(d-y)/(d-c) + b*(y-c)/(d-c);
fv = feval(f,x);
pass(1) = norm( gv - fv, inf) < f.onefun.vscale*tol

% back
f = changeMap(g, dom1);
y = ((b-a)/2) * rand(100, 1) + (b+a)/2;
fv = feval(f,y);
x = c*(b-y)/(b-a) + d*(y-a)/(b-a);
gv = feval(g,x);
pass(2) = norm( fv - gv, inf) < f.onefun.vscale*tol