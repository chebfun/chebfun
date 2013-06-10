% Test file for bndfun/flipud.m

function pass = test_flipud(pref)

if ( nargin < 1 )
    pref = bndfun.pref;
end
pref = chebtech.pref(pref);
tol = pref.bndfun.eps;
dom = [-2 7];
pass = zeros(1, 4); % Pre-allocate pass matrix

% Try some standard calls to flipud
f = bndfun(@(x) sin(x), dom, [], [], pref);
g = bndfun(@(x) sin(5-x), dom, [], [], pref);
h = flipud(f);
pass(1) = norm(g.onefun.values - h.onefun.values, inf) < f.onefun.vscale*tol;

f = bndfun(@(x) [sin(x), exp(x)], dom, [], [], pref);
g = bndfun(@(x) [sin(5-x), exp(5-x)], dom, [], [], pref);
h = flipud(f);
pass(2) = norm(g.onefun.values - h.onefun.values, inf) < 2*max(f.onefun.vscale)*tol;

f = bndfun(@(x) sin(1i*x+.5), dom, [], [], pref);
g = bndfun(@(x) sin(1i*(5-x)+.5), dom, [], [], pref);
h = flipud(f);
pass(3) = norm(g.onefun.values - h.onefun.values, inf) < 2*max(f.onefun.vscale)*tol;

f = bndfun(@(x) [sin(x), exp(1i*x)], dom, [], [], pref);
g = bndfun(@(x) [sin(5-x), exp(1i*(5-x))], dom, [], [], pref);
h = flipud(f);
pass(4) = norm(g.onefun.values - h.onefun.values, inf) < 2*max(f.onefun.vscale)*tol;

end
