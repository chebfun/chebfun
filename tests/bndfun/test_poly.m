% Test file for bndfun/poly.m

function pass = test_poly(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref;
end

% Set the domain
dom = [-2 7];

pass = zeros(1, 5); % Pre-allocate pass matrix

%%
% Check a few simple examples.

f = bndfun(@(x) zeros(size(x)), dom, [], [], pref);
p = poly(f);
pass(1) = (norm(p, inf) < f.onefun.vscale*f.onefun.epslevel);

f = bndfun(@(x) 3*ones(size(x)), dom, [], [], pref);
p = poly(f);
pass(2) = (norm(p - 3, inf) < f.onefun.vscale*f.onefun.epslevel);

f = bndfun(@(x) 6.4*x - 3i, dom, [], [], pref);
p = poly(f);
pass(3) = (norm(p - [6.4 (-3i)], inf) < f.onefun.vscale*f.onefun.epslevel);

f = bndfun(@(x) 2i*x.^5 - 3.2*x.^4 + 2*x.^2 - (1.2 + 3i), dom, [], [], pref);
p = poly(f);
pass(4) = (norm(p - [2i (-3.2) 0 2 0 -(1.2 + 3i)], inf) ...
    < f.onefun.vscale*f.onefun.epslevel);

%%
% Verify operation for array-valued bndfun objects.

f = bndfun(@(x) [3*ones(size(x)), (6.4*x - 3i), ...
    (4*x.^2 - 2i*x + 3.7)], dom, [], [], pref);
p = poly(f);
p_exact = [0 0     3;
    0 6.4   (-3i);
    4 (-2i) 3.7];
pass(5) = (norm(p(:) - p_exact(:), inf) < max(f.onefun.vscale)*f.onefun.epslevel);

end