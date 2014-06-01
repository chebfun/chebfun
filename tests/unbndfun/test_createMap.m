% Test file for unbndfun/createMap.m.

function pass = test_createMap(pref)

tol = 1e1*eps;

% Doubly-unbounded domain:
dom = [-Inf Inf];
map = unbndfun.createMap(dom);

pass(1) = ( all(map.For([-1 1]) == [-Inf Inf]) ) && ...
    ( all(map.Inv([-1e100 1e100]) - [-1 1] < tol) );

% Singly-unbounded domain (left):
dom = [-Inf 3];
map = unbndfun.createMap(dom);
pass(2) = ( all(map.For([-1 1]) == [-Inf 3]) ) && ...
    ( all(map.Inv([-1e100 3]) == [-1 1]) );

% Singly-unbounded domain (right):
dom = [-100 Inf];
map = unbndfun.createMap(dom);
pass(3) = ( all(map.For([-1 1]) == [-100 Inf]) ) && ...
    ( all(map.Inv([-100 1e100]) == [-1 1]) );

end
