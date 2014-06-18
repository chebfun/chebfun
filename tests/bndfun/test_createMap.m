% Test file for bndfun/createMap.m.

function pass = test_createMap(pref)

dom = [-2 7];
map = bndfun.createMap(dom);

pass(1) = all(map.For([-1 1]) == dom) && all(map.Inv(dom) == [-1 1]);

pass(2) = map.Der(1) == diff(dom)/2;

end
